function [t_path,x1_path,x2_path,x3_path,ttau_,iota_] = sim_cdiff_3D(x0,i0,x1_,x2_,x3_,t_,b,s,DELTA,seed)

    rng(seed)

    x1_min = x1_(1);     x1_max = x1_(end);
    x2_min = x2_(1);     x2_max = x2_(end);
    x3_min = x3_(1);     x3_max = x3_(end);

    T = t_(end);             % terminal horizon for simulation (same as the solver)
    N = round(T/DELTA);      % number of simulation steps

    t_path = linspace(0,T,N+1);             % time buckets for simulation

    x1_init = x0(1);
    x1_path = [x1_init zeros(1,N)];         % x1 path ( to be simulated )
    x2_init = x0(2);
    x2_path = [x2_init zeros(1,N)];         % x2 path ( to be simulated )
    x3_init = x0(3);
    x3_path = [x3_init zeros(1,N)];         % x3 path ( to be simulated )

    % initialize control sequences
    ttau_ = [];              % switching times
    iota_ = [];              % regime switches

    i = i0;

    steps_per_day = round( 1 / (t_(2)-t_(1)) );

    prev_day = 0;

    for n = 1:N

        x1_curr = x1_path(n);
        x2_curr = x2_path(n);
        x3_curr = x3_path(n);

        t_curr = t_path(n);

        curr_day = floor((n-1)/steps_per_day) + 1;

        if curr_day ~= prev_day
            filename = fullfile(pwd,'policy_table',sprintf('day_%03d.mat',curr_day));
            tmp = load(filename);
            ptab_day = tmp.ptab_day;
            clear tmp
        end

        inday_idx = n - (curr_day-1)*steps_per_day;        % position of the current iterate n inside the current day
                  
                                                                                    % extract best action from current policy table
        action = get_action_3D(ptab_day(:,:,:,:,inday_idx),i,x1_curr,x2_curr,x3_curr,x1_,x2_,x3_);
        if action ~= 0                                                              % if switch ...
            i = action;                                                             % ... update current regime
            ttau_ = horzcat(ttau_,t_curr);                                          % ... update control
            iota_ = horzcat(iota_,i);                                               % "                "
        end

        % simulation
        [drift1,drift2,drift3] = b{i}(t_curr,x1_curr,x2_curr,x3_curr);
        [s11,s12,s13, ...
         s21,s22,s23, ...
         s31,s32,s33] = s{i}(t_curr,x1_curr,x2_curr,x3_curr);

        dW = randn(3,1)*sqrt(DELTA);
        dx = [drift1; drift2; drift3] * DELTA + ...
             [s11 s12 s13;
              s21 s22 s23;
              s31 s32 s33] * dW;

        x1_next = x1_curr + dx(1);
        x1_next = min(max(x1_next, x1_min), x1_max);           % project x1 back onto the space grid 
        x1_path(n+1) = x1_next;
        x2_next = x2_curr + dx(2);
        x2_next = min(max(x2_next, x2_min), x2_max);           % project x2 back onto the space grid 
        x2_path(n+1) = x2_next;
        x3_next = x3_curr + dx(3);
        x3_next = min(max(x3_next, x3_min), x3_max);           % project x3 back onto the space grid 
        x3_path(n+1) = x3_next;

        prev_day = curr_day;

    end

end
