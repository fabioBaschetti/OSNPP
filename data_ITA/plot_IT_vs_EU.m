clear
close all
clc


delta = 15/(60*24);


P = [1/4 1/3 1/2 1 7/2 7 365/4 365/2 365];       % periods (in days)
omega = 2*pi ./ P;                               % frequencies

kappa_IT = 0.350010322216108;
sigma_IT = 0.111417572966839;
b_IT     = 0.611841394373854;
c_IT     = [0.410052426945638	0.160558469292021	-2.42377036825523	-1.51010901660680	0.084112231159954	0.2984233222537860	-0.0113108198750741	 0.0562621276399995	  0.0912004672700916];
s_IT     = [0.271422536765792  -0.640072655650156	 2.81556280010702	-0.95222108696697  -0.247940550661697	0.0982448206348856	-0.0161971505489276	 0.0450803116146585	 -0.0527038711819384];

phi_IT   = exp(-kappa_IT*delta);
v_IT     = sigma_IT^2/(2*kappa_IT) * (1 - phi_IT^2);


n_countries = 5;

rho      = 0;

kappa_EU = kappa_IT;
sigma_EU = sigma_IT * sqrt(n_countries*(1+(n_countries-1)*rho));
b_EU     = b_IT * (n_countries - 3);
c_EU     = c_IT * n_countries;
s_EU     = s_IT * n_countries;

phi_EU   = exp(-kappa_EU*delta);
v_EU     = sigma_EU^2/(2*kappa_EU) * (1 - phi_EU^2);


T  = 365;                              % days
dt = 1/(24*4);
Nt = T/dt;                          
t_ = linspace(0,T,Nt+1);


shift = 7;         % ... days       

t_IT = t_';
t_EU = mod(t_-shift,365)';


csin_wt_IT = cos(t_IT * omega);    
ssin_wt_IT = sin(t_IT * omega);

P_j_IT = ( (kappa_IT*csin_wt_IT(2:Nt,:) + ones(Nt-1,1)*omega.*ssin_wt_IT(2:Nt,:)) - phi_IT*(kappa_IT*csin_wt_IT(1:Nt-1,:) + ones(Nt-1,1)*omega.*ssin_wt_IT(1:Nt-1,:)) ) ./ ( ones(Nt-1,1)*(kappa_IT^2 + omega.^2) );
Q_j_IT = ( (kappa_IT*ssin_wt_IT(2:Nt,:) - ones(Nt-1,1)*omega.*csin_wt_IT(2:Nt,:)) - phi_IT*(kappa_IT*ssin_wt_IT(1:Nt-1,:) - ones(Nt-1,1)*omega.*csin_wt_IT(1:Nt-1,:)) ) ./ ( ones(Nt-1,1)*(kappa_IT^2 + omega.^2) );


csin_wt_EU = cos(t_EU * omega);    
ssin_wt_EU = sin(t_EU * omega);

P_j_EU = ( (kappa_EU*csin_wt_EU(2:Nt,:) + ones(Nt-1,1)*omega.*ssin_wt_EU(2:Nt,:)) - phi_EU*(kappa_EU*csin_wt_EU(1:Nt-1,:) + ones(Nt-1,1)*omega.*ssin_wt_EU(1:Nt-1,:)) ) ./ ( ones(Nt-1,1)*(kappa_EU^2 + omega.^2) );
Q_j_EU = ( (kappa_EU*ssin_wt_EU(2:Nt,:) - ones(Nt-1,1)*omega.*csin_wt_EU(2:Nt,:)) - phi_EU*(kappa_EU*ssin_wt_EU(1:Nt-1,:) - ones(Nt-1,1)*omega.*csin_wt_EU(1:Nt-1,:)) ) ./ ( ones(Nt-1,1)*(kappa_EU^2 + omega.^2) );


x0_IT = 0.50;
x0_EU = 0.50 * (n_countries - 3);

x_IT = zeros(1,Nt);
x_EU = zeros(1,Nt);

x_IT(1) = x0_IT;
x_EU(1) = x0_EU;
for i = 2:Nt
    mu_IT = phi_IT*x_IT(i-1) + (1-phi_IT)*b_IT + kappa_IT*(P_j_IT(i-1,:)*c_IT' + Q_j_IT(i-1,:)*s_IT');
    W_IT = randn;
    x_IT(i) = mu_IT + sqrt(v_IT)*W_IT;
    mu_EU = phi_EU*x_EU(i-1) + (1-phi_EU)*b_EU + kappa_EU*(P_j_EU(i-1,:)*c_EU' + Q_j_EU(i-1,:)*s_EU');
    W_EU = rho*W_IT + sqrt(1-rho^2)*randn;
    x_EU(i) = mu_EU + sqrt(v_EU)*W_EU;
end

% figure
% yline(0,'Color','black','LineWidth',1.5)
% hold on
% yline(1,'Color','black','LineWidth',1.5)
% hold on
% yline(n_countries-2,'Color','black','LineWidth',1.5)
% hold on
% plot((1:N)*delta,x_EU)
% hold on
% plot((1:N)*delta,x_IT)
% grid on
% xlim([0,365])

P_max = 0.9;
a = P_max*(n_countries-2);

figure
hold on

h1 = plot((1:Nt)*delta, x_EU, ...
    'Color','black', 'LineWidth',0.5);

h2 = plot((1:Nt)*delta, x_IT, ...
    'Color','red', 'LineWidth',0.5);

yline(a,     '--b', 'LineWidth',0.5);
yline(P_max, '--b', 'LineWidth',0.5);
yline(0,     '--b', 'LineWidth',0.5);

grid on

ax = gca;

% Store the limits of the left y-axis
yyaxis left
yl = ylim;

% Configure the right y-axis using the same limits
yyaxis right
ylim(yl)

% Tick values must be in increasing order
tickValues = [0, P_max, a];
tickLabels = {'$0$', '$P_{\max}$', '$a$'};

[tickValues, idx] = sort(tickValues);
tickLabels = tickLabels(idx);

yticks(tickValues)
yticklabels(tickLabels)

ax.YAxis(2).Color = 'b';
ax.YAxis(2).TickLabelInterpreter = 'latex';

% Return to the left axis for any subsequent plotting
yyaxis left

xlim([0,365])

legend([h1 h2], {'$M_t$', '$Y_t$'}, ...
    'Interpreter','latex', ...
    'Location','southeast')
