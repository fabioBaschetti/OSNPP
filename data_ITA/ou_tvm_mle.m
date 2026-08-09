function out = ou_tvm_mle(x,delta,omega)

    x = x(:);                          % as a column vector
    N = numel(x);                      % 
    n = N-1;                           % = numel(y) --> number of observations 
    t  = (0:n)'*delta;

    omega = omega(:)';                 % as a row    vector
    J = numel(omega);                  % number of ...

    csin_wt = cos(t * omega);          % (n+1) x J
    ssin_wt = sin(t * omega);          % (n+1) x J      


    % ----- profile search over \phi
    obj = @(phi) sse_at_phi(phi,delta,omega,csin_wt,ssin_wt,x);

    % tau_m = delta/3;                    
    % kappa_M = 1/tau_m;                 
    % phi_m = exp(-kappa_M*delta);
    % clear tau_m kappa_M
    % 
    % tau_M = 365/4;                     
    % kappa_m = 1/tau_M;                 
    % phi_M = exp(-kappa_m*delta);
    % clear tau_M kappa_m

    phi_m = 0+1e-05;
    phi_M = 1-1e-05;
    
    %       proceeds by identifying the optimal \phi^\star  
    [phi, ~] = fminbnd(@(ppp) obj(ppp), phi_m, phi_M);
    % -----

    %       and (conditional) \beta 
    [sse, cache] = sse_at_phi(phi,delta,omega,csin_wt,ssin_wt,x);
    
    kappa = cache.kappa;
    beta  = cache.beta;
    b = beta(1);
    c = beta(2:(J+1));
    s = beta(J+2:end);
    %       so as to minimize the sum of squared residuals
    %       eps = y - Z\beta (see AUX)

    
    % ----- MLE for the innovation variance v 
    v_hat = sse/n;
    %       and back to \sigma
    sigma = sqrt((2*kappa) / (1 - phi^2) * v_hat); 
    % -----


    out = struct();
    out.kappa = kappa; 
    out.sigma = sigma; 
    out.b = b; 
    out.c = c; 
    out.s = s;


    Z     = cache.Z;


    % ----- (naive) standard errors
    % se(\kappa)
    l_prof = @(s) -0.5*n * log(s/n);             % s = sse_\phi
    h = 1e-06*phi;
    [s_plus, ~] = sse_at_phi(phi+h,delta,omega,csin_wt,ssin_wt,x);
    [s_mins, ~] = sse_at_phi(phi-h,delta,omega,csin_wt,ssin_wt,x); 
    d2ldphi = (l_prof(s_plus)-2*l_prof(sse)+l_prof(s_mins)) / h^2;
    var_phi = -1 / d2ldphi;
    se_kappa = sqrt(var_phi / (phi*delta)^2);

    % se(\beta  | \kappa)   
    ZtZ = Z.'*Z;     R = chol(ZtZ);
    inv_ZtZ = R\(R'\eye(size(R)));
    se_beta = sqrt(diag(v_hat*inv_ZtZ));

    % se(\sigma | \kappa)
    var_sigma2 = (2*kappa/(1-phi^2))^2 * 2*v_hat^2/n;
    se_sigma = sqrt(var_sigma2) / (2*sigma);
    % -----

    out.se_kappa = se_kappa;
    out.se_sigma = se_sigma;
    out.se_b = se_beta(1);
    out.se_c = se_beta(2:(J+1));
    out.se_s = se_beta(J+2:end);

end


function [sse,cache] = sse_at_phi(phi,delta,omega,csin_wt,ssin_wt,x)

    % CONDITIONAL ON \phi (i.e. \kappa)... y_i := x_i - \phi x_{i-1} is
    % amenable to treatment by OLS regression in \beta = [b, c_j, s_j]

    kappa = -log(phi) / delta;

    N = numel(x);
    y = x(2:N) - phi*x(1:N-1);

    P_j = ( (kappa*csin_wt(2:N,:) + ones(N-1,1)*omega.*ssin_wt(2:N,:)) - phi*(kappa*csin_wt(1:N-1,:) + ones(N-1,1)*omega.*ssin_wt(1:N-1,:)) ) ./ ( ones(N-1,1)*(kappa^2 + omega.^2) );
    Q_j = ( (kappa*ssin_wt(2:N,:) - ones(N-1,1)*omega.*csin_wt(2:N,:)) - phi*(kappa*ssin_wt(1:N-1,:) - ones(N-1,1)*omega.*csin_wt(1:N-1,:)) ) ./ ( ones(N-1,1)*(kappa^2 + omega.^2) );

    Z = [(1-phi)*ones(N-1,1), kappa*P_j, kappa*Q_j];

    [Q,R] = qr(Z,0);
    beta_hat = R \ (Q'*y);   % rmk. OLS = MLE for a linear Gaussian model
    res = y - Z*beta_hat;
    sse = sum(res.^2);

    cache.kappa = kappa; 
    cache.beta  = beta_hat; 
    cache.Z = Z; 
    
end
