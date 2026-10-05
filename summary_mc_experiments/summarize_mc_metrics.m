function row = summarize_mc_metrics(label,metrics)
%SUMMARIZE_MC_METRICS One-row table of means and Monte Carlo standard errors.

    names={'TotalCost','RunningCost','SwitchingCost','AvgAbsTrackingError', ...
           'ShortageEnergy','ExcessEnergy','ShortageTimeSharePct','ExcessTimeSharePct', ...
           'NumSwitches','TimeAtBoundsPct'};
    row=table(string(label),'VariableNames',{'Policy'});
    for k=1:numel(names)
        x=metrics.(names{k});
        row.([names{k} '_Mean'])=mean(x,'omitnan');
        row.([names{k} '_SE'])=std(x,0,'omitnan')/sqrt(sum(isfinite(x)));
    end
end
