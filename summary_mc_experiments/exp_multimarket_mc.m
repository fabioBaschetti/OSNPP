function out = exp_multimarket_mc(root_dir,cfg_closed,p_closed,out_dir,varargin)
%EXP_MULTIMARKET_MC Multi-path counterpart of the paper's one-week examples.
%
% The same local residual-demand paths are used in all three environments;
% the price-taker and price-maker also use the same market residual-demand
% paths. The experiment is intended to support qualitative comparisons of
% tracking, switching and transaction volumes. Total cost levels should NOT
% be interpreted as a welfare ranking across environments because the running
% cost functions differ across the three economic environments.

    ip=inputParser;
    addParameter(ip,'Npaths',250,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'Seed',20261201,@(x)isnumeric(x)&&isscalar(x));
    parse(ip,varargin{:});

    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    cfg_pt=build_open_model(root_dir,'pricetaker');
    cfg_pm=build_open_model(root_dir,'pricemaker');
    [Y,M,diag]=simulate_open_exogenous_mc(cfg_pt,ip.Results.Npaths,ip.Results.Seed);

    % Closed economy: same Y paths, clipped to the closed-economy Y grid.
    Yc=min(max(Y,cfg_closed.x2_(1)),cfg_closed.x2_(end));
    sim_c=simulate_optimal_closed_mc(cfg_closed,p_closed,Yc);
    met_c=evaluate_closed_mc(cfg_closed,sim_c,Yc);
    row_c=summarize_mc_metrics('Closed',met_c);
    row_c.Properties.VariableNames{1}='Scenario';

    % Add open-market-only columns as NaN so rows concatenate cleanly.
    openOnly={'PurchasesEnergy','SalesEnergy','PurchaseCost','SalesRevenue','AverageSellPrice'};
    for k=1:numel(openOnly)
        row_c.([openOnly{k} '_Mean'])=NaN;
        row_c.([openOnly{k} '_SE'])=NaN;
    end

    sim_pt=simulate_optimal_open_mc(cfg_pt,Y,M);
    met_pt=evaluate_open_mc(cfg_pt,sim_pt,Y,M);
    row_pt=summarize_open_metrics('Price taker',met_pt);

    sim_pm=simulate_optimal_open_mc(cfg_pm,Y,M);
    met_pm=evaluate_open_mc(cfg_pm,sim_pm,Y,M);
    row_pm=summarize_open_metrics('Price maker',met_pm);

    % Ensure identical column ordering before vertical concatenation.
    row_c=row_c(:,row_pt.Properties.VariableNames);
    summary=[row_c;row_pt;row_pm];
    writetable(summary,fullfile(out_dir,'multimarket_mc.csv'));

    out=struct('Summary',summary,'ExogenousDiagnostics',diag, ...
               'ClosedMetrics',met_c,'PriceTakerMetrics',met_pt,'PriceMakerMetrics',met_pm);
    save(fullfile(out_dir,'multimarket_mc.mat'),'out','-v7.3');

    fprintf('\nMulti-market Monte Carlo summary\n');
    disp(summary(:,{'Scenario','AvgAbsTrackingError_Mean','NumSwitches_Mean', ...
                    'ShortageEnergy_Mean','ExcessEnergy_Mean'}));
end
