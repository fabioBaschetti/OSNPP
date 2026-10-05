function out = exp_policy_benchmark_closed(cfg,p,out_dir,varargin)
%EXP_POLICY_BENCHMARK_CLOSED Out-of-sample policy comparison in closed economy.
%
% The optimal switching policy is compared with:
%   (i) hold output constant at the common initial level;
%   (ii) a ramp-limited deadband tracking rule.
%
% The deadband is tuned on an independent training sample, then all policies
% are evaluated on the same out-of-sample residual-demand paths.

    ip=inputParser;
    addParameter(ip,'Ntrain',150,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'Ntest',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'TrainSeed',20260930,@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'TestSeed',20261001,@(x)isnumeric(x)&&isscalar(x));
    addParameter(ip,'Deadbands',0:0.01:0.20,@isnumeric);
    parse(ip,varargin{:});

    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    [Ytrain,diag_train]=simulate_residual_demand_mc(cfg,ip.Results.Ntrain,ip.Results.TrainSeed);
    hs=ip.Results.Deadbands(:);
    mean_cost=zeros(numel(hs),1);
    se_cost=zeros(numel(hs),1);
    mean_switches=zeros(numel(hs),1);

    for q=1:numel(hs)
        sim=simulate_deadband_closed_mc(cfg,Ytrain,hs(q));
        met=evaluate_closed_mc(cfg,sim,Ytrain);
        mean_cost(q)=mean(met.TotalCost);
        se_cost(q)=std(met.TotalCost)/sqrt(numel(met.TotalCost));
        mean_switches(q)=mean(met.NumSwitches);
    end

    [~,best_idx]=min(mean_cost);
    best_h=hs(best_idx);
    tuning=table(hs,mean_cost,se_cost,mean_switches, ...
        'VariableNames',{'Deadband','MeanTotalCost','SETotalCost','MeanNumSwitches'});
    writetable(tuning,fullfile(out_dir,'deadband_tuning_closed.csv'));

    [Ytest,diag_test]=simulate_residual_demand_mc(cfg,ip.Results.Ntest,ip.Results.TestSeed);

    sim_opt=simulate_optimal_closed_mc(cfg,p,Ytest);
    met_opt=evaluate_closed_mc(cfg,sim_opt,Ytest);

    sim_db=simulate_deadband_closed_mc(cfg,Ytest,best_h);
    met_db=evaluate_closed_mc(cfg,sim_db,Ytest);

    sim_hold=simulate_hold_closed_mc(cfg,Ytest);
    met_hold=evaluate_closed_mc(cfg,sim_hold,Ytest);

    summary=[summarize_mc_metrics('Optimal switching',met_opt); ...
             summarize_mc_metrics(sprintf('Deadband h=%.3f',best_h),met_db); ...
             summarize_mc_metrics('Hold P=P_0',met_hold)];

    vopt=summary.TotalCost_Mean(1);
    summary.CostGapVsOptimalPct=100*(summary.TotalCost_Mean-vopt)/vopt;

    writetable(summary,fullfile(out_dir,'policy_benchmark_closed.csv'));

    out=struct();
    out.BestDeadband=best_h;
    out.Tuning=tuning;
    out.Summary=summary;
    out.TrainDiagnostics=diag_train;
    out.TestDiagnostics=diag_test;
    out.MetricsOptimal=met_opt;
    out.MetricsDeadband=met_db;
    out.MetricsHold=met_hold;
    save(fullfile(out_dir,'policy_benchmark_closed.mat'),'out','-v7.3');

    fprintf('\nClosed-economy policy benchmark\n');
    fprintf('Best training deadband: %.4f\n',best_h);
    disp(summary(:,{'Policy','TotalCost_Mean','TotalCost_SE','CostGapVsOptimalPct', ...
                    'AvgAbsTrackingError_Mean','NumSwitches_Mean'}));
end
