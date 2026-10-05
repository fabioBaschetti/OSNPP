function out = exp_mesh_stability_closed(root_dir,out_dir)
%EXP_MESH_STABILITY_CLOSED Coarse/baseline/fine numerical stability check.
%
% The physical ramping speed is fixed across meshes. The experiment reports
% value-function stability at representative states and policy disagreement
% on a common evaluation lattice at t=0 and t=T/2.

    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    addpath(fullfile(root_dir,'utils'));
    addpath(fullfile(root_dir,'utils','2D'));

    variants=struct( ...
        'Name', {'Coarse','Baseline','Fine'}, ...
        'Nx1', {81,101,121}, ...
        'Nx2', {109,136,163}, ...
        'dt', {1/72,1/96,1/120});

    states=[0.35 0.45; 0.50 0.50; 0.75 0.65];
    evalP=linspace(0.20,0.90,41);
    evalY=linspace(-0.10,1.15,51);

    nV=numel(variants);
    vals=zeros(nV,size(states,1),3);
    maps0=cell(nV,1);
    mapsMid=cell(nV,1);
    runSeconds=zeros(nV,1);

    for q=1:nV
        fprintf('\nMesh stability: solving %s grid...\n',variants(q).Name);
        cfg=build_closed_model('Nx1',variants(q).Nx1,'Nx2',variants(q).Nx2, ...
                               'dt',variants(q).dt,'T',7);
        tic;
        [v,p]=SL2D_OS_FH(cfg.x1_,cfg.x2_,cfg.t_,cfg.b,cfg.s,cfg.f,cfg.h,cfg.g,'mode','min');
        runSeconds(q)=toc;

        for s=1:size(states,1)
            for i=1:3
                vals(q,s,i)=interp2(cfg.x1_,cfg.x2_,v(:,:,i,1),states(s,1),states(s,2),'linear');
            end
        end

        maps0{q}=extract_policy_map(p,cfg,evalP,evalY,0);
        mapsMid{q}=extract_policy_map(p,cfg,evalP,evalY,cfg.T/2);
        clear v p cfg
    end

    ref=vals(end,:,:);
    maxAbs=zeros(nV,1); meanAbs=zeros(nV,1);
    dis0=zeros(nV,1); disMid=zeros(nV,1);
    for q=1:nV
        d=abs(vals(q,:,:)-ref);
        maxAbs(q)=max(d(:));
        meanAbs(q)=mean(d(:));
        dis0(q)=100*mean(maps0{q}(:)~=maps0{end}(:));
        disMid(q)=100*mean(mapsMid{q}(:)~=mapsMid{end}(:));
    end

    summary=table(string({variants.Name})',[variants.Nx1]',[variants.Nx2]',[variants.dt]', ...
        runSeconds,maxAbs,meanAbs,dis0,disMid, ...
        'VariableNames',{'Mesh','NxP','NxY','dtDays','RuntimeSeconds','MaxAbsValueDiffVsFine', ...
                         'MeanAbsValueDiffVsFine','PolicyDisagreement_t0_Pct','PolicyDisagreement_tmid_Pct'});

    meshCol=strings(nV*size(states,1)*3,1);
    pCol=zeros(size(meshCol)); yCol=pCol; regimeCol=pCol; valueCol=pCol;
    row=0;
    for q=1:nV
        for s=1:size(states,1)
            for i=1:3
                row=row+1;
                meshCol(row)=string(variants(q).Name);
                pCol(row)=states(s,1); yCol(row)=states(s,2);
                regimeCol(row)=i; valueCol(row)=vals(q,s,i);
            end
        end
    end
    valuesLong=table(meshCol,pCol,yCol,regimeCol,valueCol, ...
        'VariableNames',{'Mesh','P','Y','RegimeIndex','Value'});

    writetable(summary,fullfile(out_dir,'mesh_stability_closed.csv'));
    writetable(valuesLong,fullfile(out_dir,'mesh_values_closed.csv'));
    out=struct('Summary',summary,'Values',valuesLong,'Variants',variants,'States',states);
    save(fullfile(out_dir,'mesh_stability_closed.mat'),'out','-v7.3');

    fprintf('\nMesh stability summary\n');
    disp(summary);
end

function A=extract_policy_map(p,cfg,evalP,evalY,t)
    [PP,YY]=meshgrid(evalP,evalY);
    A=zeros(numel(evalY),numel(evalP),3,'uint8');
    tidx=max(1,min(cfg.Nt,round(t/cfg.dt)+1));
    k1=nearest_uniform_index(PP,cfg.x1_);
    k2=nearest_uniform_index(YY,cfg.x2_);
    for i=1:3
        ii=i*ones(size(PP));
        tt=tidx*ones(size(PP));
        idx=sub2ind(size(p),k2,k1,ii,tt);
        raw=uint8(p(idx));
        raw(raw==0)=uint8(i); % compare post-decision regime, not 0 coding
        A(:,:,i)=raw;
    end
end

function idx=nearest_uniform_index(y,xgrid)
    dx=(xgrid(end)-xgrid(1))/(numel(xgrid)-1);
    idx=round(1+(y-xgrid(1))/dx);
    idx=max(1,min(numel(xgrid),idx));
end
