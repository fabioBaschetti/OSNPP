clear
close all
clc


PARENT_DIR = fileparts(pwd);

UTL_gn_PATH = PARENT_DIR  + "\" + "utils";
UTL_2D_PATH = UTL_gn_PATH + "\" + "2D";
addpath(UTL_gn_PATH)
addpath(UTL_2D_PATH)


%% -------------------------------- set-up --------------------------------

P_min    = 0.20;
P_max    = 0.90;

% --- space grids
x1_min = P_min;
x1_max = P_max;
Nx1 = 101;  
x1_ = linspace(x1_min, x1_max, Nx1);   % x1 = P
x2_min = -0.15;
x2_max =  1.20;
Nx2 = 136;
x2_ = linspace(x2_min, x2_max, Nx2);   % x2 = \bar{Y}


% --- time grid
T  =   7;                              % days
dt = 1/(24*4);
Nt = T/dt;                          
t_ = linspace(0,T,Nt+1);


% --- regimes
r =  0.05 / dt;                         % ramping: 5% / 15-minutes      
Ivals   = [-r, 0, +r];

m = numel(Ivals);


% --- coefficients of the controlled dynamics

P = [1/4 1/3 1/2 1 7/2 7 365/4 365/2 365];       % periods (in days)
omega = 2*pi ./ P;                               % frequencies

% residual demand (IT)
kappa = 0.350010322216108;
mu    = 0.611841394373854;
xi    = [0.410052426945638	0.160558469292021	-2.42377036825523	-1.51010901660680	0.084112231159954	0.2984233222537860	-0.0113108198750741	 0.0562621276399995	  0.0912004672700916];
zeta  = [0.271422536765792  -0.640072655650156	 2.81556280010702	-0.95222108696697  -0.247940550661697	0.0982448206348856	-0.0161971505489276	 0.0450803116146585	 -0.0527038711819384];
sigma = 0.111417572966839;


% drift
theta = @(t) mu + sum(xi.*cos(omega*t) + zeta.*sin(omega*t));
b = cell(m,1); 
for k = 1:m
    b{k} = @(t,x,y) deal(max(Ivals(k),0)*(x<P_max) + min(Ivals(k),0)*(x>P_min), kappa*(theta(t)-y)); 
end
clear k

% diffu 
s = cell(m,1);
for k = 1:m
    s{k} = @(t,x,y) deal(0*x,         0*x, ...
                         0*y, sigma + 0*y);
end
clear k


% --- running costs
delta       = 0.08;

                                       % price when the (EU-)clearing source is
s_res = 0.00;                          % renewables
s_npp = 0.20;                          % nuclear
s_fss = 0.40;                          % fossil fuels

gamma1  = s_npp + 0.50*delta;
lambda2 = (s_fss + delta);
lambda1 = s_res;

f = cell(1,m);
for k = 1:m
    f{k} = @(t,x,y) lambda1*max(x-y,0) + lambda2*max(y-x,0) + gamma1*x;
end

% --- terminal cost 
h = cell(1,m);
for k = 1:m
    h{k} = @(x,y) 0.00*x + 0.00*y;
end


% --- switching costs
G = zeros(m);
G(1,2) = 0.0020; G(2,1) = 0.0008;
G(2,3) = 0.0024; G(3,2) = 0.0002;
G(1,3) = 0.0035; G(3,1) = 0.0008;
G = 0.20 * G;
g = @(i,j,t,x,y) G(i,j) + 0.00*t + 0.00*x + 0.00*y;


%% -------------------------- numerical solution --------------------------

[v,p] = SL2D_OS_FH(x1_,x2_,t_,b,s,f,h,g,'mode','min');

t_idx = 1;
regime = ["-r", "0", "+r"];
for k = 1:m
    
    figure
    M = p(:,:,k,t_idx);
    M(M==0) = k;
    h = imagesc(x1_,x2_,M);                 
    set(h,'Interpolation','nearest')
    colormap([0 0 1;         % -r -> blue
              0 1 0;         %  0 -> green
              1 0 0]);       % +r -> red
    axis image
    xlabel('$P$', 'Interpreter','latex')
    ylabel('$Y$', 'Interpreter','latex')
    title(sprintf('$I(t=0,(P,Y),i_{0^-}=%s)$', regime(k)), 'Interpreter','latex')

    hold on

    ax = gca;
    xl = xlim(ax);
    yl = ylim(ax);
    
    P_point = 0.5;
    Y_point = 0.5;
    plot(P_point, Y_point, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6)
    plot([xl(1) P_point], [Y_point Y_point], 'k:', 'LineWidth',1.5)
    plot([P_point P_point], [Y_point yl(2)], 'k:', 'LineWidth',1.5)

end
clear k


%% ------------------------------ simulation ------------------------------ 

seed = 35;

x1_init = 0.5;
x2_init = 0.5;

x0 = [x1_init,x2_init];
i0 = 2;

DELTA = dt;


do_sim = false;

if do_sim

    [t_path,x1_path,x2_path,ttau_,iota_] = sim_cdiff_2D(x0,i0,x1_,x2_,t_,p,b,s,DELTA,seed);

else

    TRJ_PATH = PARENT_DIR  + "\" + "paths";
    addpath(TRJ_PATH)

    load(sprintf("x2_path_%d.mat", seed))   % x2_path

    [t_path,x1_path,ttau_,iota_] = sim_conditional_2D(x0,i0,x1_,x2_,t_,p,b,DELTA,x2_path);

end


%% -------------------------------- plots ---------------------------------

YP_clsdmkt_plt(t_path,x1_path,x2_path,P_min,P_max)

if isempty(ttau_)
    disp('no switches')
else
    switches_plt(i0,T,ttau_,iota_,m)
end
