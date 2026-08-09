clear
close all
clc


PARENT_DIR = fileparts(pwd);

UTL_gn_PATH = PARENT_DIR  + "\" + "utils";
UTL_3D_PATH = UTL_gn_PATH + "\" + "3D";
addpath(UTL_gn_PATH)
addpath(UTL_3D_PATH)


%% -------------------------------- set-up --------------------------------

P_min_IT = 0.20;
P_max_IT = 0.90;

% --- space grids
x1_min = P_min_IT;
x1_max = P_max_IT;
Nx1 = 76; 
x1_ = linspace(x1_min, x1_max, Nx1);   % x1 = P^{IT}
x2_min = -0.20;
x2_max =  1.40;
Nx2 = 161;
x2_ = linspace(x2_min, x2_max, Nx2);   % x2 = Y^{IT}

n = 5;                                 % countries in the EU (including IT)

x3_min = -1.5;
x3_max =  4.0;
Nx3 = 138;
x3_ = linspace(x3_min, x3_max, Nx3);   % x3 = Y^{EU}


% --- time grid
T  =   7;                              % days
dt = 1/(24*4);
Nt = T/dt;                          
t_ = linspace(0,T,Nt+1);


% --- regimes
r = 0.05 / dt;                         % ramping: 5% / 15-minutes      
Ivals   = [-r, 0, +r];

m = numel(Ivals);


% --- coefficients of the controlled dynamics

P = [1/4 1/3 1/2 1 7/2 7 365/4 365/2 365];       % periods (in days)
omega = 2*pi ./ P;                               % frequencies

% residual demand (IT)
kappa_IT = 0.350010322216108;
mu_IT    = 0.611841394373854;
xi_IT    = [0.410052426945638	0.160558469292021	-2.42377036825523	-1.51010901660680	0.084112231159954	0.2984233222537860	-0.0113108198750741	 0.0562621276399995	  0.0912004672700916];
zeta_IT  = [0.271422536765792  -0.640072655650156	 2.81556280010702	-0.95222108696697  -0.247940550661697	0.0982448206348856	-0.0161971505489276	 0.0450803116146585	 -0.0527038711819384];
sigma_IT = 0.111417572966839;


% residual demand (EU)
kappa_EU = kappa_IT;
mu_EU    = mu_IT * (n*0.40);
xi_EU    = xi_IT * n;
zeta_EU  = zeta_IT * n;
rho      = 0;                          % correlation between the Brownian shocks in EU and IT
sigma_EU = sigma_IT*sqrt(n*(1+(n-1)*rho));


% drift
b_params.Ivals    = Ivals;
b_params.P_min_IT = P_min_IT;
b_params.P_max_IT = P_max_IT;
b_params.kappa_IT = kappa_IT;
b_params.kappa_EU = kappa_EU;
theta_IT = zeros(1,Nt);
theta_EU = zeros(1,Nt);
% delay the EU process by "shift" ...
shift = 7;         % ... days       
for k = 1:Nt
    tk_IT = t_(k);
    theta_IT(k) = mu_IT + sum(xi_IT.*cos(omega*tk_IT) + zeta_IT.*sin(omega*tk_IT));
    tk_EU = mod(t_(k)-shift,365);
    theta_EU(k) = mu_EU + sum(xi_EU.*cos(omega*tk_EU) + zeta_EU.*sin(omega*tk_EU));
end
b_params.theta_IT = theta_IT;
b_params.theta_EU = theta_EU;

% diffu 
C = zeros(3,3);
C(2,2) =               sigma_IT;
C(3,2) =           rho*sigma_EU;
C(3,3) = sqrt(1-rho^2)*sigma_EU;


% --- running costs
delta       = 0.08;                    % penalization for buying from the market

                                       % price when the (EU-)clearing source is
s_res = 0.00;                          % renewables
s_npp = 0.20;                          % nuclear
s_fss = 0.40;                          % fossil fuels

P_max_EU = (n*0.60)*P_max_IT;                 


gamma1  = s_npp + 0.50*delta;
f = @(x,y,z) -  phi(z+y-x,P_max_EU,s_res,s_npp,s_fss)       .*max(x-y,0) ...
             + (phi(z+y-x,P_max_EU,s_res,s_npp,s_fss)+delta).*max(y-x,0) ...
             + gamma1*x;


% --- switching costs
G = zeros(m);
G(1,2) = 0.0020; G(2,1) = 0.0008;
G(2,3) = 0.0024; G(3,2) = 0.0002;
G(1,3) = 0.0035; G(3,1) = 0.0008;
G = 0.20 * G;


%% -------------------------- numerical solution --------------------------

do_run = false;

if do_run
    v0 = SL3D_OS_FH(x1_,x2_,x3_,t_,b_params,C,f,G,'min');
else
    % do nothing
end


%% ------------------------------ simulation ------------------------------

seed = 35;

x1_init = 0.5;
x2_init = 0.5;
x3_init = x2_init*(n*0.40);

x0 = [x1_init,x2_init,x3_init];
i0 = 2;

DELTA = dt;

b = cell(m,1);
for k = 1:m
    b{k} = @(t,x,y,z) deal(max(Ivals(k),0)*(x<P_max_IT) + min(Ivals(k),0)*(x>P_min_IT), ...
                           kappa_IT*(mu_IT+sum(xi_IT.*cos(omega*t               )+zeta_IT.*sin(omega*t               ))-y), ...
                           kappa_EU*(mu_EU+sum(xi_EU.*cos(omega*mod(t-shift,365))+zeta_EU.*sin(omega*mod(t-shift,365)))-z) );
end
clear k

s = cell(m,1);
for k = 1:m
    s{k} = @(t,x,y,z) deal(0*x,          0*x,          0*x, ...
                           0*y, C(2,2) + 0*y,          0*y, ...
                           0*z, C(3,2) + 0*z, C(3,3) + 0*z);
end
clear k


do_sim = false;

if do_sim

    [t_path,x1_path,x2_path,x3_path,ttau_,iota_] = sim_cdiff_3D(x0,i0,x1_,x2_,x3_,t_,b,s,DELTA,seed);

else
    
    TRJ_PATH = PARENT_DIR  + "\" + "paths";
    addpath(TRJ_PATH)

    load(sprintf("x2_path_%d.mat", seed))   % x2_path
    load(sprintf("x3_path_%d.mat", seed))   % x3_path

    [t_path,x1_path,ttau_,iota_] = sim_conditional_3D(x0,i0,x1_,x2_,x3_,t_,b,DELTA,x2_path,x3_path);

end

z_path = x3_path + x2_path - x1_path;

S_path = nan(1,Nt+1);
S_path(z_path <= 0) = s_res;
S_path(z_path > P_max_EU) = s_fss;
S_path(isnan(S_path)) = s_npp;


%% -------------------------------- plots --------------------------------- 

price_determinant_plt(P_max_EU,t_path,z_path,'$Z_t$')

YP_openmkt_plt(t_path,x1_path,x2_path,P_min_IT,P_max_IT,P_max_EU,z_path)

if isempty(ttau_)
    disp('no switches')
else
    switches_plt(i0,T,ttau_,iota_,m)
end


% ------------------------------ AUXILIARIES ------------------------------
function out = phi(x,threshold,y_low,y_med,y_hig)    
    out = y_med*ones(size(x));
    out(x <= 0) = y_low;
    out(x > threshold) = y_hig;
end
% -------------------------------------------------------------------------
