function sim = simulate_optimal_open_mc(cfg,Y,M)
%SIMULATE_OPTIMAL_OPEN_MC Apply stored 3D policy tables to common paths.
%
% Policy files are loaded one day at a time so that the large 3D policy
% tables do not all need to be resident in memory simultaneously.

    Npaths=size(Y,1);
    Nt=cfg.Nt;
    if size(Y,2)~=Nt+1 || size(M,2)~=Nt+1
        error('Y and M must have Nt+1 columns.');
    end
    if ~exist(cfg.PolicyDir,'dir')
        error('Policy directory not found: %s',cfg.PolicyDir);
    end

    P=zeros(Npaths,Nt+1);
    P(:,1)=cfg.x1_init;
    regime=cfg.i0*ones(Npaths,1);
    num_switches=zeros(Npaths,1);
    switch_cost=zeros(Npaths,1);
    steps_per_day=round(1/cfg.dt);

    current_day=-1;
    ptab_day=[];

    for n=1:Nt
        day=floor((n-1)/steps_per_day)+1;
        if day~=current_day
            fname=fullfile(cfg.PolicyDir,sprintf('day_%03d.mat',day));
            if ~exist(fname,'file')
                error(['Missing policy file %s. If the repository was cloned with Git LFS, ' ...
                       'run git lfs pull before this experiment.'],fname);
            end
            tmp=load(fname,'ptab_day');
            if ~isfield(tmp,'ptab_day')
                error('File %s does not contain ptab_day (possible unresolved Git-LFS pointer).',fname);
            end
            ptab_day=tmp.ptab_day;
            current_day=day;
        end

        inday=n-(day-1)*steps_per_day;
        slice=ptab_day(:,:,:,:,inday);
        k1=nearest_uniform_index(P(:,n),cfg.x1_);
        k2=nearest_uniform_index(Y(:,n),cfg.x2_);
        k3=nearest_uniform_index(M(:,n),cfg.x3_);
        idx=sub2ind(size(slice),k2,k1,k3,regime);
        action=double(slice(idx));

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
