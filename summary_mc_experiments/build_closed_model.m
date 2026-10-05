function cfg = build_closed_model(varargin)
%BUILD_CLOSED_MODEL Parameters for the closed-economy benchmark.
%
% cfg = build_closed_model('Nx1',101,'Nx2',136,'dt',1/96,'T',7)
%
% The physical ramping speed is fixed at 4.8 normalized power units/day,
% corresponding to 5 percentage points over 15 minutes. It is deliberately
% kept independent of the numerical time step so that mesh-refinement tests
% do not change the underlying operating model.

    p = inputParser;
    addParameter(p,'Nx1',101,@(x)isnumeric(x)&&isscalar(x)&&x>=3);
    addParameter(p,'Nx2',136,@(x)isnumeric(x)&&isscalar(x)&&x>=3);
    addParameter(p,'dt',1/(24*4),@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p,'T',7,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(p,varargin{:});

    cfg.P_min = 0.20;
    cfg.P_max = 0.90;

    cfg.x1_ = linspace(cfg.P_min,cfg.P_max,p.Results.Nx1);
    cfg.x2_ = linspace(-0.15,1.20,p.Results.Nx2);

    cfg.T  = p.Results.T;
    cfg.dt = p.Results.dt;
    cfg.Nt = round(cfg.T/cfg.dt);
    cfg.t_ = linspace(0,cfg.T,cfg.Nt+1);
    cfg.dt = cfg.t_(2)-cfg.t_(1); % exact mesh actually used

    % Physical ramping parameter: 5 percentage points in 15 minutes.
    cfg.r = 0.05/(1/(24*4));
    cfg.Ivals = [-cfg.r,0,+cfg.r];
    cfg.m = numel(cfg.Ivals);

    periods = [1/4 1/3 1/2 1 7/2 7 365/4 365/2 365];
    cfg.omega = 2*pi./periods;

    cfg.kappa = 0.350010322216108;
    cfg.mu     = 0.611841394373854;
    cfg.xi     = [0.410052426945638 0.160558469292021 -2.42377036825523 -1.51010901660680 ...
                  0.084112231159954 0.2984233222537860 -0.0113108198750741 0.0562621276399995 0.0912004672700916];
    cfg.zeta   = [0.271422536765792 -0.640072655650156 2.81556280010702 -0.95222108696697 ...
                 -0.247940550661697 0.0982448206348856 -0.0161971505489276 0.0450803116146585 -0.0527038711819384];
    cfg.sigma  = 0.111417572966839;

    cfg.theta = @(t) cfg.mu + sum(cfg.xi.*cos(cfg.omega*t) + cfg.zeta.*sin(cfg.omega*t));

    cfg.b = cell(cfg.m,1);
    for k=1:cfg.m
        ik = cfg.Ivals(k);
        cfg.b{k} = @(t,x,y) deal(max(ik,0).*(x<cfg.P_max) + min(ik,0).*(x>cfg.P_min), ...
                                  cfg.kappa.*(cfg.theta(t)-y));
    end

    cfg.s = cell(cfg.m,1);
    for k=1:cfg.m
        cfg.s{k} = @(t,x,y) deal(0*x,0*x,0*y,cfg.sigma+0*y);
    end

    cfg.delta = 0.08;
    cfg.s_res = 0.00;
    cfg.s_npp = 0.20;
    cfg.s_fss = 0.40;

    cfg.gamma1  = cfg.s_npp + 0.50*cfg.delta;
    cfg.lambda2 = cfg.s_fss + cfg.delta;
    cfg.lambda1 = cfg.s_res;

    cfg.f = cell(1,cfg.m);
    for k=1:cfg.m
        cfg.f{k} = @(t,x,y) cfg.lambda1.*max(x-y,0) + cfg.lambda2.*max(y-x,0) + cfg.gamma1.*x;
    end

    cfg.h = cell(1,cfg.m);
    for k=1:cfg.m
        cfg.h{k} = @(x,y) 0*x + 0*y;
    end

    G = zeros(cfg.m);
    G(1,2)=0.0020; G(2,1)=0.0008;
    G(2,3)=0.0024; G(3,2)=0.0002;
    G(1,3)=0.0035; G(3,1)=0.0008;
    cfg.G = 0.20*G;
    cfg.g = @(i,j,t,x,y) cfg.G(i,j) + 0*t + 0*x + 0*y;

    cfg.x1_init = 0.5;
    cfg.x2_init = 0.5;
    cfg.i0 = 2;
end
