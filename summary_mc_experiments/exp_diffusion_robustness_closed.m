function out = exp_diffusion_robustness_closed(cfg,p,out_dir,deadband,varargin)
%EXP_DIFFUSION_ROBUSTNESS_CLOSED Stress-test the baseline policy.
%
% The policy is NOT re-optimized under the alternative data-generating
% processes. This is an out-of-sample parameter-misspecification exercise.
% Every scenario remains a Brownian Ornstein--Uhlenbeck diffusion of the
% same form as in the paper. We vary only the volatility and mean-reversion
% speed of the residual-demand process. The same Gaussian shock matrix is
% used across scenarios (common random numbers) so that scenario differences
% are not confounded by different Monte Carlo draws.

    if nargin<4 || isempty(deadband), deadband=0.05; end

    ip=inputParser;
    addParameter(ip,'Npaths',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'BaseSeed',20261001,@(x)isnumeric(x)&&isscalar(x));
    parse(ip,varargin{:});

    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    sc=struct( ...
        'Name', { ...
            'Lower volatility: 0.75 sigma', ...
            'Baseline diffusion', ...
            'Higher volatility: 1.25 sigma', ...
            'Slower mean reversion: 0.75 kappa', ...
            'Faster mean reversion: 1.25 kappa'}, ...
        'SigmaMult', {0.75, 1.00, 1.25, 1.00, 1.00}, ...
        'KappaMult', {1.00, 1.00, 1.00, 0.75, 1.25});

    rows=table();
    all_metrics=cell(numel(sc),2);
    diagnostics=cell(numel(sc),1);

    for q=1:numel(sc)
        [Y,diag]=simulate_residual_demand_mc(cfg,ip.Results.Npaths,ip.Results.BaseSeed, ...
            'sigma_multiplier',sc(q).SigmaMult, ...
            'kappa_multiplier',sc(q).KappaMult);
        diagnostics{q}=diag;

        sim_opt=simulate_optimal_closed_mc(cfg,p,Y);
        met_opt=evaluate_closed_mc(cfg,sim_opt,Y);
        r1=summarize_mc_metrics('Baseline optimal policy',met_opt);
        r1.Scenario=string(sc(q).Name);
        r1.SigmaMultiplier=sc(q).SigmaMult;
        r1.KappaMultiplier=sc(q).KappaMult;
        r1.YGridBoundaryHitPct=diag.boundary_hit_share_pct;

        sim_db=simulate_deadband_closed_mc(cfg,Y,deadband);
        met_db=evaluate_closed_mc(cfg,sim_db,Y);
        r2=summarize_mc_metrics(sprintf('Deadband h=%.3f',deadband),met_db);
        r2.Scenario=string(sc(q).Name);
        r2.SigmaMultiplier=sc(q).SigmaMult;
        r2.KappaMultiplier=sc(q).KappaMult;
        r2.YGridBoundaryHitPct=diag.boundary_hit_share_pct;

        all_metrics{q,1}=met_opt;
        all_metrics{q,2}=met_db;
        rows=[rows; r1; r2]; %#ok<AGROW>
    end

    rows=movevars(rows,{'Scenario','SigmaMultiplier','KappaMultiplier'},'Before','Policy');
    writetable(rows,fullfile(out_dir,'diffusion_robustness_closed.csv'));

    out=struct('Summary',rows,'Diagnostics',{diagnostics},'Metrics',{all_metrics},'Deadband',deadband);
    save(fullfile(out_dir,'diffusion_robustness_closed.mat'),'out','-v7.3');

    fprintf('\nClosed-economy diffusion-parameter robustness\n');
    disp(rows(:,{'Scenario','Policy','TotalCost_Mean','TotalCost_SE', ...
                 'AvgAbsTrackingError_Mean','NumSwitches_Mean','YGridBoundaryHitPct'}));
end
