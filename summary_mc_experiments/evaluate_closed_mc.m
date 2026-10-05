function metrics = evaluate_closed_mc(cfg,sim,Y)
%EVALUATE_CLOSED_MC Pathwise economic and operating metrics.

    P=sim.P;
    dt=cfg.dt;
    P0=P(:,1:end-1);
    Y0=Y(:,1:end-1);
    shortage=max(Y0-P0,0);
    excess=max(P0-Y0,0);

    running=dt*sum(cfg.lambda1*excess + cfg.lambda2*shortage + cfg.gamma1*P0,2);
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
end
