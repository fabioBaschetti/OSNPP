function v0 = SL3D_OS_FH(x1_,x2_,x3_,t_,b_params,C,f,G,mode)

    policy_dir = fullfile(pwd,'policy_table');

    if exist(policy_dir,'dir')
        rmdir(policy_dir,'s');
    end
    mkdir(policy_dir);

    switch lower(mode)

        case 'max'
            v0 = solver(x1_,x2_,x3_,t_,b_params,C,f,G,policy_dir);

        case 'min'
            % flip sign of ...
            % ... running  rewards
            f_hat = f;
            f = @(x,y,z) -f_hat(x,y,z);
            clear f_hat
            % solve max problem
            v0 = solver(x1_,x2_,x3_,t_,b_params,C,f,G,policy_dir);
            % mirror value function
            v0 = -v0;

        otherwise
            error('ERROR: invalid mode');

    end

end


function v0 = solver(x1_,x2_,x3_,t_,b_params,C,f,G,policy_dir)

    % --- get time-grid information
    Nt = length(t_) - 1;
    steps_per_day = round( 1 / (t_(2)-t_(1)) );
    % ---


    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool;
    end
    n_cores = pool.NumWorkers;


    Nx1 = length(x1_);
    Nx2 = length(x2_);
    Nx3 = length(x3_);


    n_chunks = max(1, min(n_cores, round(Nx3/35)));
    chunks   = make_chunks(x3_, n_chunks);

    
    % --- prepare chunks of the state variables
    X1_chunks = cell(1,n_chunks);
    X2_chunks = cell(1,n_chunks);
    X3_chunks = cell(1,n_chunks);
    parfor q = 1:n_chunks
        x3_loc = chunks{q}.x3_;
        [X1_loc,X2_loc,X3_loc] = meshgrid(x1_,x2_,x3_loc);
        X1_chunks{q} = X1_loc;
        X2_chunks{q} = X2_loc;
        X3_chunks{q} = X3_loc;
    end
    % ---


    % --- precompute running rewards on chunks
    %     (this is only possible because f does not depend on time)
    f_chunks = cell(1,n_chunks);
    parfor q = 1:n_chunks
        X1_loc = X1_chunks{q};
        X2_loc = X2_chunks{q};
        X3_loc = X3_chunks{q};  
        f_chunks{q} = f(X1_loc,X2_loc,X3_loc);
    end
    % ---

    
    % --- store grid information for linear interpolation
    grid.x1_offset = x1_(1);
    grid.x2_offset = x2_(1);
    grid.x3_offset = x3_(1);

    grid.Nx1 = Nx1;
    grid.Nx2 = Nx2;
    grid.Nx3 = Nx3;
                  
                                            % h = x_(2) - x(1) is ...
    grid.inv_h1 = 1 / (x1_(2) - x1_(1));    % (uniform) step size on the x1_in grid
    grid.inv_h2 = 1 / (x2_(2) - x2_(1));    % (uniform) step size on the x2_in grid
    grid.inv_h3 = 1 / (x3_(2) - x3_(1));    % (uniform) step size on the x3_in grid
    % ---


    m  = length(b_params.Ivals);       % number of regimes [ I = 1 ,..., m ]


    % --- extract drift parameters  
    Ivals    = b_params.Ivals;
    P_min_IT = b_params.P_min_IT;
    P_max_IT = b_params.P_max_IT;
    kappa_IT = b_params.kappa_IT;
    kappa_EU = b_params.kappa_EU;
    theta_IT = b_params.theta_IT;
    theta_EU = b_params.theta_EU;
    % ---


    % --- pre-compute \sum_{i=1}^d e_1 \sigma_{i1} + \dots + e_d \sigma_{id}
    %     to be used for the footpoints
    %     N.B. this is only possible because the diffusion matrix is constant
    %          (it neither depends on time nor on the space variables)

    % diffusion coefficients
    s11 = C(1,1);  s12 = C(1,2);  s13 = C(1,3);
    s21 = C(2,1);  s22 = C(2,2);  s23 = C(2,3);
    s31 = C(3,1);  s32 = C(3,2);  s33 = C(3,3);

    % sign matrix for the footpoints   
    % n_row = 2^3
    % n_col = 3   [e1,e2,e3]
    sgn = dec2bin(0:7) - '0';
    sgn = 2*sgn - 1;

    sig = zeros(8,3);
    for k = 1:8
        
        e1 = sgn(k,1);
        e2 = sgn(k,2);
        e3 = sgn(k,3);

        sig(k,1) = e1*s11 + e2*s12 + e3*s13;          
        sig(k,2) = e1*s21 + e2*s22 + e3*s23;
        sig(k,3) = e1*s31 + e2*s32 + e3*s33;

    end

    % ---


    % --- Floyd–Warshall algorithm for all-pairs shortest path matrix D
    D = G;                                  % D is constant through time --> it can be pre-computed outside the time loop 
    for l = 1:m                             % and over the spatial nodes --> we just need the mxm matrix (regime pairs) 
        D = min(D, D(:,l) + D(l,:));        %                                rather than a d+2-dimensional tensor with           
    end                                     %                                the whole spacial lattice
    % ---


    % ---------- terminal condition ---------- 
    v_next = zeros(Nx2,Nx1,Nx3,m);          % v^i_T = 0 in the absence of terminal costs (h = 0)
    

    prev_day = floor((Nt-1)/steps_per_day) + 1;
    % --- initialize policy table for last day
    ptab_day = zeros(Nx2,Nx1,Nx3,m,steps_per_day,'uint8');
    % ---


    % ---------- backward recursion ----------
    for n = Nt:-1:1

        t_curr = t_(n);
        dt = t_(n+1) - t_curr;

        sqdt = sqrt(dt);

        thetaIT_n = theta_IT(n);
        thetaEU_n = theta_EU(n);

        v_pieces = cell(1,n_chunks);
        p_pieces = cell(1,n_chunks);

        parfor q = 1:n_chunks

            idx_inchunk = chunks{q}.idx;
            d_chunk = length(idx_inchunk);

            X1_loc = X1_chunks{q};
            X2_loc = X2_chunks{q};
            X3_loc = X3_chunks{q};

            % --- build continuation candidates C_i^{(n)}(x_,v_i^{(n+1)})
            %                                                __.__.__.__
            %     N.B. v_i^{(n+1)} is known from previous step at t + dt
            continuation = zeros(Nx2,Nx1,d_chunk,m);

            f_loc = f_chunks{q};

            % regime-independent drift components
            b2_tmp = kappa_IT * (thetaIT_n - X2_loc);
            b3_tmp = kappa_EU * (thetaEU_n - X3_loc);

            for i = 1:m

                % regime-dependent drift component
                b1_tmp = max(Ivals(i),0)*(X1_loc<P_max_IT) + min(Ivals(i),0)*(X1_loc>P_min_IT);

                cum = zeros(Nx2,Nx1,d_chunk);
                for k = 1:8

                    Y1 = X1_loc + b1_tmp*dt + sig(k,1)*sqdt;
                    Y2 = X2_loc + b2_tmp*dt + sig(k,2)*sqdt;
                    Y3 = X3_loc + b3_tmp*dt + sig(k,3)*sqdt;

                    cum = cum + my_interp3D(v_next(:,:,:,i),grid,Y1,Y2,Y3);

                end

                continuation(:,:,:,i) = f_loc*dt + 0.125*cum;

            end
            % ---


            % --- per-slice solution for the value function via shortest paths

            % apply path formula: v_i(.) = max_j (C_j(.) - D_ij(.))
            v_loc = -Inf(Nx2,Nx1,d_chunk,m);
            for i = 1:m
                for j = 1:m
                    candidate = continuation(:,:,:,j) - D(i,j);
                    v_loc(:,:,:,i) = max(v_loc(:,:,:,i), candidate);
                end
            end

            % ---


            % --- extract policy table
            %     according to
            %     deterministic tie rules:
            %          - "continue" vs " switch " ties [C   == S  ] --> prioritize continuation 
            %          - " switch " vs " switch " ties [S_i == S_j] --> prioritize i < j  (*)
            %            (*) assumption: lower regimes are more 'economic' and to be preferred

            atol = 1e-10;
            rtol = 1e-10;

            p_loc = zeros(Nx2,Nx1,d_chunk,m,'uint8');

            for i = 1:m
                optimal_switch_val = -Inf*ones(Nx2,Nx1,d_chunk);               % initialize the switching operations at -Inf (we take the max! (°))
                optimal_switch_reg = zeros(Nx2,Nx1,d_chunk);
                for j = 1:m
                    if j == i
                        continue
                    end
                    cndidat_switch_val = v_loc(:,:,:,j) - G(i,j);                                  % j_th candidate for switching
                    target = optimal_switch_val;
                    target(~isfinite(target))=0;
                    scl_ss = max(abs(cndidat_switch_val),abs(target));
                    tol_ss = atol + rtol*scl_ss;
                    mask = cndidat_switch_val > (optimal_switch_val+tol_ss);                       %                                           (°)                                              
                                                                                                   % if the current switch is convenient relative to the previous ones ...
                    optimal_switch_val(mask)  = cndidat_switch_val(mask);                          % ... update the value
                    optimal_switch_reg(mask)  = j;                                                 % ... update the regime
                end
                scl_cs = max(abs(continuation(:,:,:,i)),abs(optimal_switch_val));
                tol_cs = atol + rtol*scl_cs;
                aux = optimal_switch_reg;                                              % compare best switch with continuation
                aux(continuation(:,:,:,i)>=(optimal_switch_val-tol_cs)) = 0;
                p_loc(:,:,:,i) = uint8(aux);
            end

            % ---

            v_pieces{q} = v_loc;
            p_pieces{q} = p_loc;

        end

        v_curr = zeros(Nx2,Nx1,Nx3,m);
        p_curr = zeros(Nx2,Nx1,Nx3,m,'uint8');

        for q = 1:n_chunks
            idx_inchunk = chunks{q}.idx;
            v_curr(:,:,idx_inchunk,:) = v_pieces{q};
            p_curr(:,:,idx_inchunk,:) = p_pieces{q};
        end

        % --- write policy table day by day
        
        curr_day = floor((n-1)/steps_per_day) + 1;
        
        if curr_day ~= prev_day
            
            filename = fullfile(policy_dir,sprintf('day_%03d.mat',prev_day));
            save(filename,'ptab_day','-v7.3');

            ptab_day = zeros(Nx2,Nx1,Nx3,m,steps_per_day,'uint8');
            prev_day = curr_day;

        end

        n_low = (curr_day-1)*steps_per_day + 1;
        inday_idx = n - n_low + 1;          % position of the current iterate n inside the current day

        ptab_day(:,:,:,:,inday_idx) = p_curr; 

        % ---

        v_next = v_curr;

    end

    % --- write policy table for day 1
    filename = fullfile(policy_dir,sprintf('day_%03d.mat',prev_day));
    save(filename,'ptab_day','-v7.3');
    % ---

    v0 = v_next;

    % ----------------------------------------

end


function chunks = make_chunks(x3_,n_chunks)

    Nx3 = length(x3_);

    buckts = round(linspace(0,Nx3,n_chunks+1));
    chunks = cell(1,n_chunks);

    for q = 1:n_chunks
        l = buckts(q) + 1;
        h = buckts(q+1);
        chunks{q}.idx = l:h;
        chunks{q}.x3_ = x3_(l:h);
    end

end
