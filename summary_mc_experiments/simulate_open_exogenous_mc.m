function [Y,M,diag] = simulate_open_exogenous_mc(cfg,Npaths,seed)
%SIMULATE_OPEN_EXOGENOUS_MC Common local and market residual-demand paths.
%
% The discretization follows the repository's 3D simulator. The same draws
% can be used for both price-taker and price-maker policy evaluation.

    rng(seed,'twister');
    Nt=cfg.Nt; dt=cfg.dt;
    Y=zeros(Npaths,Nt+1); M=zeros(Npaths,Nt+1);
    Y(:,1)=cfg.x2_init; M(:,1)=cfg.x3_init;
    hitY=0; hitM=0;

    for n=1:Nt
        t=cfg.t_(n);
        e1=randn(Npaths,1);
        e2=randn(Npaths,1);
        dWY=e1;
        dWM=cfg.rho*e1 + sqrt(1-cfg.rho^2)*e2;

        ynew=Y(:,n)+cfg.kappa_IT*(cfg.theta_IT(t)-Y(:,n))*dt + cfg.sigma_IT*sqrt(dt).*dWY;
        mnew=M(:,n)+cfg.kappa_EU*(cfg.theta_EU(t)-M(:,n))*dt + cfg.sigma_EU*sqrt(dt).*dWM;

        hitY=hitY+sum(ynew<cfg.x2_(1) | ynew>cfg.x2_(end));
        hitM=hitM+sum(mnew<cfg.x3_(1) | mnew>cfg.x3_(end));
        Y(:,n+1)=min(max(ynew,cfg.x2_(1)),cfg.x2_(end));
        M(:,n+1)=min(max(mnew,cfg.x3_(1)),cfg.x3_(end));
    end

    diag=struct('YBoundaryHitPct',100*hitY/(Npaths*Nt), ...
                'MBoundaryHitPct',100*hitM/(Npaths*Nt), ...
                'Seed',seed,'Npaths',Npaths);
end
