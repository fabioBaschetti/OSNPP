function [v,p] = SL2D_OS_FH_evo_full(x1_,x2_,t_,b,s,f,h,g,varargin)

    % ---------- defaults ----------
    mode = 'max';
    per_slice_sln = 'shortest-path';
    extrap = 'clamp';

    for k = 1:2:length(varargin)
        if k+1 >length(varargin)
            error('Missing value for property: %s',varargin{k})
        else
            switch varargin{k}
                case 'mode'
                    mode = varargin{k+1};
                case 'per_slice_sln'
                    per_slice_sln = varargin{k+1};
                case 'extrap'
                    extrap = varargin{k+1};
                otherwise
                    error('Invalid property name: %s',varargin{k})
            end
        end
    end

    switch lower(mode)
    
        case 'max'
            [v,p] = solver(x1_,x2_,t_,b,s,f,h,g,per_slice_sln,extrap);
                           
        case 'min'
            % flip sign of ...
            % ... running  rewards
            f = cellfun(@(fun) @(varargin) -fun(varargin{:}), f, 'UniformOutput', false);
            % ... terminal rewards
            h = cellfun(@(fun) @(varargin) -fun(varargin{:}), h, 'UniformOutput', false);
            % solve max problem
            [v,p] = solver(x1_,x2_,t_,b,s,f,h,g,per_slice_sln,extrap);
            % mirror value function
            v = -v;

        otherwise
            error('ERROR: invalid mode');

    end

end

function [v,p] = solver(x1_,x2_,t_,b,s,f,h,g,per_slice_sln,extrap)

    % solve the 'max' problem

    Nx1 = length(x1_);
    Nx2 = length(x2_);
    
    Nt = length(t_) - 1;

    m = length(b);                % number of regimes [ I = 1 ,..., m ]


    % build 2-dim space grid (Nx2 x Nx1)
    [X1,X2] = meshgrid(x1_, x2_);  


    % initialize value and policy tables 
    v = zeros(Nx2,Nx1,m,Nt+1);
    p = zeros(Nx2,Nx1,m,Nt  );


    % ---------- terminal condition ---------- 
    for i = 1:m
        v(:,:,i,Nt+1) = h{i}(X1,X2);
    end
    clear i


    % ---------- backward recursion ----------
    for n = Nt:-1:1
    
        t_curr = t_(n);
        dt = t_(n+1) - t_curr;

        % build continuation candidates C_i^{(n)}(x_,v_i^{(n+1)})
        %                                            __.__.__.__
        % N.B. v_i^{(n+1)} is known from previous step at t + dt
        continuation = zeros(Nx2,Nx1,m);
        for i = 1:m

            % of x1  of x2
            [b1_tmp,b2_tmp] = b{i}(t_curr,X1,X2);                            % tmp = step  : n
            [s11_tmp,s12_tmp,s21_tmp,s22_tmp] = s{i}(t_curr,X1,X2);          %       regime: i

            % get (temporary) footpoints for linear interpolation
            % of the value function [v_i^{(n+1)}]
            y1_pp = X1 + b1_tmp*dt + (+s11_tmp+s12_tmp)*sqrt(dt);
            y1_pm = X1 + b1_tmp*dt + (+s11_tmp-s12_tmp)*sqrt(dt);
            y1_mp = X1 + b1_tmp*dt + (-s11_tmp+s12_tmp)*sqrt(dt);
            y1_mm = X1 + b1_tmp*dt + (-s11_tmp-s12_tmp)*sqrt(dt);
            y2_pp = X2 + b2_tmp*dt + (+s21_tmp+s22_tmp)*sqrt(dt);
            y2_pm = X2 + b2_tmp*dt + (+s21_tmp-s22_tmp)*sqrt(dt);
            y2_mp = X2 + b2_tmp*dt + (-s21_tmp+s22_tmp)*sqrt(dt);
            y2_mm = X2 + b2_tmp*dt + (-s21_tmp-s22_tmp)*sqrt(dt);

            continuation(:,:,i) = f{i}(t_curr,X1,X2)*dt + 0.25*(my_interp2D(v(:,:,i,n+1),x1_,x2_,y1_pp,y2_pp,extrap) + ...
                                                                my_interp2D(v(:,:,i,n+1),x1_,x2_,y1_pm,y2_pm,extrap) + ...
                                                                my_interp2D(v(:,:,i,n+1),x1_,x2_,y1_mp,y2_mp,extrap) + ...
                                                                my_interp2D(v(:,:,i,n+1),x1_,x2_,y1_mm,y2_mm,extrap));

        end

        switch lower(per_slice_sln)
        
            case 'switching-closure'
                U = continuation;        % initialize U at continuation --> convergence
                for sweeps = 1:(m-1)     % (from below) in m -1 sweeps
                    U_old = U;           % update U through the sweeps
                    % norm of update
                    res = -Inf;          % initialize at -Inf (we take the max! (*))
                    for i = 1:m
                        best = continuation(:,:,i);
                        for j = 1:m
                            if j == i
                                continue
                            end
                            i2j_candidate = U_old(:,:,j) - g(i,j,t_curr,X1,X2);
                            best = max(best,i2j_candidate);
                        end
                        delta = max(abs(best-U_old(:,:,i)), [], 'all');      % norm of update
                        U(:,:,i) = best;
                        res = max(res,delta);                                %  (*)
                    end
                    if res == 0          % if closure stabilizes earlier
                        break   
                    end
                end

            case 'shortest-path'
                % initialize D = D_{t_curr}(X1,X2,i,j) for i,j = 1,...,m
                D = zeros(Nx2,Nx1,m,m);
                for i = 1:m
                    for j = 1:m
                        if j == i         % g(i,i,t,x1,x2) = 0 for all t \in [0,T]   
                            continue      %                    for all x \in \mathbbm{R}^2
                        end               % or self-loops are convenient
                        D(:,:,i,j) = g(i,j,t_curr,X1,X2);
                    end
                end
                % Floyd–Warshall algorithm for all-pairs shortest path matrix over
                % all spatial nodes (X1,X2)
                for l = 1:m
                    D = min(D, D(:,:,:,l) + D(:,:,l,:));
                end
                % apply path formula: U_i(.) = max_j (C_j(.) - D_ij(.)) 
                C = repmat(reshape(continuation, Nx2,Nx1,1,m), [1 1 m 1]);  % 4D view of 'continuation' with a dummy 
                                                                            % i-dimension to broadcast against D
                T = C - D;                                                  % T(i,j,Nx2,Nx1) = C_j - D_ij
                U = max(T,[],4);                                            % max over j (dimension 2)
                
            otherwise
                error('ERROR: per-slice solver %s is not supported.', per_slice_sln);

        end
        v_curr = U;
        v(:,:,:,n) = v_curr;

        % extract policy table according to ...
        % ... deterministic tie rules:
        %     - "continue" vs " switch " ties [C   == S  ] --> prioritize continuation 
        %     - " switch " vs " switch " ties [S_i == S_j] --> prioritize i < j  (*)
        %       (*) assumption: lower regimes are more 'economic' and to be preferred
        atol = 1e-10;
        rtol = 1e-10;
        for i = 1:m
            optimal_switch_val = -Inf*ones(Nx2,Nx1);       % initialize the switching operations at -Inf (we take the max! (°))
            optimal_switch_reg = zeros(Nx2,Nx1);
            for j = 1:m
                if j == i
                    continue
                end
                cndidat_switch_val = v_curr(:,:,j) - g(i,j,t_curr,X1,X2);      % j_th candidate for switching
                target = optimal_switch_val;
                target(~isfinite(target))=0;
                scl_ss = max(abs(cndidat_switch_val),abs(target));
                tol_ss = atol + rtol*scl_ss;
                mask = cndidat_switch_val > (optimal_switch_val+tol_ss);       %                                           (°)                                              
                                                                               % if the current switch is convenient relative to the previous ones ...
                optimal_switch_val(mask)  = cndidat_switch_val(mask);          % ... update the value
                optimal_switch_reg(mask)  = j;                                 % ... update the regime
            end
            scl_cs = max(abs(continuation(:,:,i)),abs(optimal_switch_val));    
            tol_cs = atol + rtol*scl_cs;
            aux = optimal_switch_reg;                                          % compare best switch with continuation
            aux(continuation(:,:,i)>=(optimal_switch_val-tol_cs)) = 0;
            p(:,:,i,n) = aux;
        end
    
    end
    % ----------------------------------------

end
