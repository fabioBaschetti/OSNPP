function sim = simulate_deadband_closed_mc(cfg,Y,deadband)
%SIMULATE_DEADBAND_CLOSED_MC Ramp toward residual demand outside a deadband.
%
% Desired regime at each decision time:
%   +r if Y-P > deadband
%   -r if P-Y > deadband
%    0 otherwise.
% Direct changes of regime pay the same switching-cost matrix as the model.

    Npaths=size(Y,1);
    Nt=cfg.Nt;
    P=zeros(Npaths,Nt+1);
    P(:,1)=cfg.x1_init;
    regime=cfg.i0*ones(Npaths,1);
    num_switches=zeros(Npaths,1);
    switch_cost=zeros(Npaths,1);

    for n=1:Nt
        gap=Y(:,n)-P(:,n);
        desired=2*ones(Npaths,1); % zero-ramp regime
        desired(gap> deadband)=3;
        desired(gap<-deadband)=1;

        sw=(desired~=regime);
        if any(sw)
            gidx=sub2ind(size(cfg.G),regime(sw),desired(sw));
            switch_cost(sw)=switch_cost(sw)+cfg.G(gidx);
            regime(sw)=desired(sw);
            num_switches(sw)=num_switches(sw)+1;
        end

        P(:,n+1)=P(:,n)+reshape(cfg.Ivals(regime),[],1)*cfg.dt;
        P(:,n+1)=min(max(P(:,n+1),cfg.P_min),cfg.P_max);
    end

    sim=struct('P',P,'NumSwitches',num_switches,'SwitchingCost',switch_cost);
end
