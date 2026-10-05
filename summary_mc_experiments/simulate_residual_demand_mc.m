function [Y,diag] = simulate_residual_demand_mc(cfg,Npaths,seed,varargin)
%SIMULATE_RESIDUAL_DEMAND_MC Common residual-demand paths for policy tests.
%
% Optional parameters:
%   'sigma_multiplier'  default 1
%   'kappa_multiplier'  default 1
%
% All scenarios remain within the Brownian Ornstein--Uhlenbeck diffusion
% class used in the paper. Paths are Euler simulated and clamped to the
% same Y-grid used by the numerical policy, matching the treatment in the
% repository simulators.

    p=inputParser;
    addParameter(p,'sigma_multiplier',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p,'kappa_multiplier',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(p,varargin{:});

    rng(seed,'twister');
    Nt=cfg.Nt;
    dt=cfg.dt;
    Y=zeros(Npaths,Nt+1);
    Y(:,1)=cfg.x2_init;
    hit_count=0;

    sigma_eff=p.Results.sigma_multiplier*cfg.sigma;
    kappa_eff=p.Results.kappa_multiplier*cfg.kappa;

    for n=1:Nt
        t=cfg.t_(n);
        dW=randn(Npaths,1)*sqrt(dt);
        ynew=Y(:,n) + kappa_eff*(cfg.theta(t)-Y(:,n))*dt + sigma_eff*dW;
        hit_count = hit_count + sum(ynew<cfg.x2_(1) | ynew>cfg.x2_(end));
        Y(:,n+1)=min(max(ynew,cfg.x2_(1)),cfg.x2_(end));
    end

    diag=struct();
    diag.boundary_hit_share_pct=100*hit_count/(Npaths*Nt);
    diag.seed=seed;
    diag.Npaths=Npaths;
    diag.sigma_multiplier=p.Results.sigma_multiplier;
    diag.kappa_multiplier=p.Results.kappa_multiplier;
end
