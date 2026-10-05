function cfg = build_open_model(root_dir,market_type)
%BUILD_OPEN_MODEL Parameters shared by the two open-economy experiments.
%
% market_type: 'pricetaker' or 'pricemaker'

    if nargin<2, market_type='pricetaker'; end
    market_type=lower(string(market_type));
    if ~ismember(market_type,["pricetaker","pricemaker"])
        error('market_type must be pricetaker or pricemaker.');
    end

    cfg.MarketType=market_type;
    cfg.P_min=0.20;
    cfg.P_max=0.90;
    cfg.x1_=linspace(cfg.P_min,cfg.P_max,76);
    cfg.x2_=linspace(-0.20,1.40,161);
    cfg.x3_=linspace(-1.5,4.0,138);

    cfg.T=7;
    cfg.dt=1/(24*4);
    cfg.Nt=round(cfg.T/cfg.dt);
    cfg.t_=linspace(0,cfg.T,cfg.Nt+1);

    cfg.r=0.05/(1/(24*4));
    cfg.Ivals=[-cfg.r,0,+cfg.r];
    cfg.i0=2;
    cfg.x1_init=0.5;
    cfg.x2_init=0.5;

    cfg.n=5;
    cfg.x3_init=cfg.x2_init*(cfg.n*0.40);

    periods=[1/4 1/3 1/2 1 7/2 7 365/4 365/2 365];
    cfg.omega=2*pi./periods;

    cfg.kappa_IT=0.350010322216108;
    cfg.mu_IT=0.611841394373854;
    cfg.xi_IT=[0.410052426945638 0.160558469292021 -2.42377036825523 -1.51010901660680 ...
               0.084112231159954 0.2984233222537860 -0.0113108198750741 0.0562621276399995 0.0912004672700916];
    cfg.zeta_IT=[0.271422536765792 -0.640072655650156 2.81556280010702 -0.95222108696697 ...
                -0.247940550661697 0.0982448206348856 -0.0161971505489276 0.0450803116146585 -0.0527038711819384];
    cfg.sigma_IT=0.111417572966839;

    cfg.kappa_EU=cfg.kappa_IT;
    cfg.mu_EU=cfg.mu_IT*(cfg.n*0.40);
    cfg.xi_EU=cfg.xi_IT*cfg.n;
    cfg.zeta_EU=cfg.zeta_IT*cfg.n;
    cfg.rho=0;
    cfg.sigma_EU=cfg.sigma_IT*sqrt(cfg.n*(1+(cfg.n-1)*cfg.rho));
    cfg.shift=7;

    cfg.theta_IT=@(t) cfg.mu_IT + sum(cfg.xi_IT.*cos(cfg.omega*t) + cfg.zeta_IT.*sin(cfg.omega*t));
    cfg.theta_EU=@(t) cfg.mu_EU + sum(cfg.xi_EU.*cos(cfg.omega*mod(t-cfg.shift,365)) + ...
                                     cfg.zeta_EU.*sin(cfg.omega*mod(t-cfg.shift,365)));

    cfg.delta=0.08;
    cfg.s_res=0.00;
    cfg.s_npp=0.20;
    cfg.s_fss=0.40;
    cfg.P_max_EU=(cfg.n*0.60)*cfg.P_max;
    cfg.gamma1=cfg.s_npp+0.50*cfg.delta;

    G=zeros(3);
    G(1,2)=0.0020; G(2,1)=0.0008;
    G(2,3)=0.0024; G(3,2)=0.0002;
    G(1,3)=0.0035; G(3,1)=0.0008;
    cfg.G=0.20*G;

    switch market_type
        case "pricetaker"
            cfg.PolicyDir=fullfile(root_dir,'OSNPP_openmkt_pricetaker','policy_table');
        case "pricemaker"
            cfg.PolicyDir=fullfile(root_dir,'OSNPP_openmkt_pricemaker','policy_table');
    end
end
