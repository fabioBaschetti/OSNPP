clear 
close all
clc


load("total_load.mat")
load("renewables.mat")


tab = total_load;
tab.Load = tab.Load - renewables.Production;

clear renewables total_load


x = tab.Load / (max(tab.Load));
x = 1.25*x;

figure
plot(tab.DateTime,x, 'Color','red')
grid on
% ylabel('MW')
% title('Residual Demand (scaled)')


delta = 15/(60*24);

P = [1/4 1/3 1/2 1 7/2 7 365/4 365/2 365];       % periods (in days)
omega = 2*pi ./ P;                               % frequencies

mle = ou_tvm_mle(x,delta,omega);

kappa = mle.kappa;
sigma = mle.sigma;
b = mle.b;
c = mle.c;
s = mle.s;

phi = exp(-kappa*delta);
v   = sigma^2/(2*kappa) * (1 - phi^2);


N = numel(x);
n = N-1;

t = (0:n)' * delta;
csin_wt = cos(t * omega);    
ssin_wt = sin(t * omega);

P_j = ( (kappa*csin_wt(2:N,:) + ones(N-1,1)*omega.*ssin_wt(2:N,:)) - phi*(kappa*csin_wt(1:N-1,:) + ones(N-1,1)*omega.*ssin_wt(1:N-1,:)) ) ./ ( ones(N-1,1)*(kappa^2 + omega.^2) );
Q_j = ( (kappa*ssin_wt(2:N,:) - ones(N-1,1)*omega.*csin_wt(2:N,:)) - phi*(kappa*ssin_wt(1:N-1,:) - ones(N-1,1)*omega.*csin_wt(1:N-1,:)) ) ./ ( ones(N-1,1)*(kappa^2 + omega.^2) );

mu = phi*x(1:end-1) + (1-phi)*b + kappa*(P_j*c + Q_j*s);


x_hat = zeros(N,1);
x_hat(1) = x(1);
for i = 2:N
    x_hat(i) = mu(i-1) + sqrt(v)*randn;
end

hold on
oliveGreen = [0.60 0.70 0.35];
plot(tab.DateTime,x_hat, 'Color',oliveGreen)
legend({'$Y$','$\hat{Y}$'}, 'Interpreter','latex', 'Location','northwest')


figure
plot(tab.DateTime,x, 'Color','red')
hold on
plot(tab.DateTime,x_hat,'Color',oliveGreen)
xlim([ ...
    datetime('01-Mar-2024 00:00:01','InputFormat','dd-MMM-yyyy HH:mm:ss'), ...
    datetime('15-Mar-2024 23:59:59','InputFormat','dd-MMM-yyyy HH:mm:ss') ...
])
grid on
% ylabel('MW')
% title('Residual Demand (scaled)')
legend({'$Y$','$\hat{Y}$'}, 'Interpreter','latex', 'Location','southeast')


eps = x(2:end) - mu;    % innovations
z = eps / sqrt(v);


figure
histogram(z,'Normalization','pdf')
hold on
x_grid = linspace(min(z),max(z),100); 
carrot = [230 126 34]/255;
plot(x_grid, normpdf(x_grid, 0,1), 'Color',carrot, 'LineWidth',2)
clear x_grid
grid on
xlabel('z')
ylabel('pdf')


[h, p] = jbtest(z);
if h == 1
    result = 'failed';
else
    result = 'passed';
end
fprintf('Jarque-Bera test for normality %s with critical value %.4f\n', ...
        result, p);


maxLag = 4*24;
figure
autocorr(z, 'NumLags',maxLag)
title([])
ylabel('ACF($z$)', 'Interpreter','latex')

[h,p,~,~] = lbqtest(z, 'Lags',maxLag, 'Alpha',0.05);
if h == 1
    result = 'failed';
else
    result = 'passed';
end
fprintf('Ljung-Box test for residual autocorrelation %s with critical value %.4f\n', ...
        result, p);


figure
autocorr(z.^2, 'NumLags',maxLag)
title([])
ylabel('ACF($z^2$)', 'Interpreter','latex')
