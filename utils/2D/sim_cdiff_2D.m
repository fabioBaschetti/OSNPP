function [t_path,x1_path,x2_path,ttau_,iota_] = sim_cdiff_2D_OS_FH_full(x0,i0,x1_,x2_,t_,p,b,s,DELTA,seed)

    rng(seed)

    x1_min = x1_(1);     x1_max = x1_(end);
    x2_min = x2_(1);     x2_max = x2_(end);

    T = t_(end);             % terminal horizon for simulation (same as the solver)
    N = round(T/DELTA);      % number of simulation steps
    
    t_path = linspace(0,T,N+1);              % time buckets for simulation

    x1_init = x0(1);
    x1_path = [x1_init zeros(1,N)];          % x1 path ( to be simulated )
    x2_init = x0(2);
    x2_path = [x2_init zeros(1,N)];          % x2 path ( to be simulated )

    % initialize control sequences
    ttau_ = [];              % switching times
    iota_ = [];              % regime switches

    i = i0;

    for n = 1:N

        x1_curr = x1_path(n);
        x2_curr = x2_path(n);

        t_curr = t_path(n);
    
        [~,idx] = min(abs(t_(1:end-1) - t_curr));                    % index of the policy table to be used for time t_curr 

        action = 0;

        action = get_action(p(:,:,:,idx),i,x1_curr,x2_curr,x1_,x2_); % extract best action from current policy table
        if action ~= 0                                               % if switch ...
            i = action;                                              % ... update current regime
            ttau_ = horzcat(ttau_,t_curr);                           % ... update control
            iota_ = horzcat(iota_,i);                                % "                "
        end

        % simulation 
        [ drift1, drift2] = b{i}(t_curr,x1_curr,x2_curr);
        [s11,s12,s21,s22] = s{i}(t_curr,x1_curr,x2_curr);
            
        dW = randn(2,1)*sqrt(DELTA);
        dx = [drift1;drift2]*DELTA + [s11 s12; s21 s22]*dW;

        x1_next = x1_curr + dx(1);
        x1_next = min(max(x1_next, x1_min), x1_max);           % project x1 back onto the space grid 
        x1_path(n+1) = x1_next;
        x2_next = x2_curr + dx(2);
        x2_next = min(max(x2_next, x2_min), x2_max);           % project x2 back onto the space grid 
        x2_path(n+1) = x2_next;

    end

end


function out = get_nearest_grid_node(y,x_)
    Nx = length(x_);
    x0 = x_(1);
    dx = (x_(end) - x0) / (Nx - 1);
    out = max(1, min(Nx, round(1+(y-x0)/dx)));
end

function out = get_action(p,i,y1,y2,x1_,x2_)
    k1 = get_nearest_grid_node(y1,x1_);
    k2 = get_nearest_grid_node(y2,x2_);
    out = p(k2,k1,i);
end
