function metrics = evaluate_open_mc(cfg,sim,Y,M)
%EVALUATE_OPEN_MC Pathwise metrics for price-taker or price-maker policy.

    P=sim.P;
    dt=cfg.dt;
    P0=P(:,1:end-1); Y0=Y(:,1:end-1); M0=M(:,1:end-1);
    shortage=max(Y0-P0,0);
    excess=max(P0-Y0,0);

    switch cfg.MarketType
        case "pricetaker"
            Z=M0;
        case "pricemaker"
            Z=M0+Y0-P0;
        otherwise
            error('Unknown market type.');
    end

    S=cfg.s_npp*ones(size(Z));
    S(Z<=0)=cfg.s_res;
    S(Z>cfg.P_max_EU)=cfg.s_fss;

    purchase=dt*sum((S+cfg.delta).*shortage,2);
    sales=dt*sum(S.*excess,2);
    operating=dt*sum(cfg.gamma1*P0,2);
    running=purchase-sales+operating;
    total=running+sim.SwitchingCost;

    metrics=struct();
    metrics.TotalCost=total;
    metrics.RunningCost=running;
    metrics.SwitchingCost=sim.SwitchingCost;
    metrics.AvgAbsTrackingError=mean(abs(P0-Y0),2);
    metrics.ShortageEnergy=dt*sum(shortage,2);
    metrics.ExcessEnergy=dt*sum(excess,2);
    metrics.ShortageTimeSharePct=100*mean(Y0>P0,2);
    metrics.ExcessTimeSharePct=100*mean(P0>Y0,2);
    metrics.NumSwitches=sim.NumSwitches;
    tol=0.5*min(diff(cfg.x1_));
    metrics.TimeAtBoundsPct=100*mean(abs(P0-cfg.P_min)<=tol | abs(P0-cfg.P_max)<=tol,2);
    metrics.PurchasesEnergy=metrics.ShortageEnergy;
    metrics.SalesEnergy=metrics.ExcessEnergy;
    metrics.PurchaseCost=purchase;
    metrics.SalesRevenue=sales;
    metrics.AverageSellPrice=mean(S,2);
end
