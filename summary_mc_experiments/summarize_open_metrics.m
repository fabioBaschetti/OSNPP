function row = summarize_open_metrics(label,metrics)
%SUMMARIZE_OPEN_METRICS Means and Monte Carlo standard errors for open market.
    names={'TotalCost','RunningCost','SwitchingCost','AvgAbsTrackingError', ...
           'ShortageEnergy','ExcessEnergy','ShortageTimeSharePct','ExcessTimeSharePct', ...
           'NumSwitches','TimeAtBoundsPct','PurchasesEnergy','SalesEnergy', ...
           'PurchaseCost','SalesRevenue','AverageSellPrice'};
    row=table(string(label),'VariableNames',{'Scenario'});
    for k=1:numel(names)
        x=metrics.(names{k});
        row.([names{k} '_Mean'])=mean(x,'omitnan');
        row.([names{k} '_SE'])=std(x,0,'omitnan')/sqrt(sum(isfinite(x)));
    end
end
