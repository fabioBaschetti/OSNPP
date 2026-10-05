function sim = simulate_optimal_closed_mc(cfg,p,Y)
%SIMULATE_OPTIMAL_CLOSED_MC Apply the computed 2D switching policy to paths Y.

    Npaths=size(Y,1);
    Nt=cfg.Nt;
    if size(Y,2)~=Nt+1
        error('Y must have Nt+1 columns.');
    end

    P=zeros(Npaths,Nt+1);
    P(:,1)=cfg.x1_init;
    regime=cfg.i0*ones(Npaths,1);
    num_switches=zeros(Npaths,1);
    switch_cost=zeros(Npaths,1);

    for n=1:Nt
        k1=nearest_uniform_index(P(:,n),cfg.x1_);
        k2=nearest_uniform_index(Y(:,n),cfg.x2_);
        tn=n*ones(Npaths,1);
        idx=sub2ind(size(p),k2,k1,regime,tn);
        action=double(p(idx));

        sw=(action~=0);
        if any(sw)
            gidx=sub2ind(size(cfg.G),regime(sw),action(sw));
            switch_cost(sw)=switch_cost(sw)+cfg.G(gidx);
            regime(sw)=action(sw);
            num_switches(sw)=num_switches(sw)+1;
        end

        P(:,n+1)=P(:,n)+reshape(cfg.Ivals(regime),[],1)*cfg.dt;
        P(:,n+1)=min(max(P(:,n+1),cfg.P_min),cfg.P_max);
    end

    sim=struct('P',P,'NumSwitches',num_switches,'SwitchingCost',switch_cost);
end

function idx=nearest_uniform_index(y,xgrid)
    dx=(xgrid(end)-xgrid(1))/(numel(xgrid)-1);
    idx=round(1+(y-xgrid(1))/dx);
    idx=max(1,min(numel(xgrid),idx));
end
