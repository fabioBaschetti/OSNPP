function outputs = run_mc_experiments(varargin)

% Optional name/value arguments:
%   'RunPolicyBenchmark'   true/false (default true)
%   'RunRobustness'        true/false (default true)
%   'RunMultiMarketMC'     true/false (default true)
%   'RunMeshStability'     true/false (default true)
%   'Ntrain'               default 150
%   'Ntest'                default 1000
%   'Nmulti'               default 250

    ip=inputParser;
    addParameter(ip,'RunPolicyBenchmark',true,@islogical);
    addParameter(ip,'RunRobustness',true,@islogical);
    addParameter(ip,'RunMultiMarketMC',true,@islogical);
    addParameter(ip,'RunMeshStability',true,@islogical);
    addParameter(ip,'Ntrain',150,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'Ntest',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'Nmulti',250,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    parse(ip,varargin{:});

    this_dir=fileparts(mfilename('fullpath'));
    root_dir=fileparts(this_dir);
    out_dir=fullfile(this_dir,'output');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    addpath(this_dir);
    addpath(fullfile(root_dir,'utils'));
    addpath(fullfile(root_dir,'utils','2D'));
    addpath(fullfile(root_dir,'utils','3D'));

    outputs=struct();

    need_baseline = ip.Results.RunPolicyBenchmark || ip.Results.RunRobustness || ip.Results.RunMultiMarketMC;
    if need_baseline
        fprintf('\nBuilding baseline closed-economy policy once for reuse...\n');
        cfg=build_closed_model();
        [~,p_closed]=SL2D_OS_FH(cfg.x1_,cfg.x2_,cfg.t_,cfg.b,cfg.s,cfg.f,cfg.h,cfg.g,'mode','min');
    else
        cfg=[]; p_closed=[];
    end

    best_h=0.05;
    if ip.Results.RunPolicyBenchmark
        outputs.PolicyBenchmark=exp_policy_benchmark_closed(cfg,p_closed,out_dir, ...
            'Ntrain',ip.Results.Ntrain,'Ntest',ip.Results.Ntest);
        best_h=outputs.PolicyBenchmark.BestDeadband;
    end

    if ip.Results.RunRobustness
        outputs.Robustness=exp_diffusion_robustness_closed(cfg,p_closed,out_dir,best_h, ...
            'Npaths',ip.Results.Ntest);
    end

    if ip.Results.RunMultiMarketMC
        try
            outputs.MultiMarketMC=exp_multimarket_mc(root_dir,cfg,p_closed,out_dir, ...
                'Npaths',ip.Results.Nmulti);
        catch ME
            warning('Multi-market Monte Carlo skipped/failed: %s',ME.message);
            outputs.MultiMarketMCError=ME;
        end
    end

    clear p_closed

    if ip.Results.RunMeshStability
        outputs.MeshStability=exp_mesh_stability_closed(root_dir,out_dir);
    end

    save(fullfile(out_dir,'summary_mc_experiments_workspace.mat'),'outputs','-v7.3');
    fprintf('\nAll requested MC experiments completed. Output folder:\n%s\n',out_dir);
    
end
