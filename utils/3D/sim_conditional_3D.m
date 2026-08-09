function [t_path,x1_path,ttau_,iota_] = sim_conditional_3D(x0,i0,x1_,x2_,x3_,t_,b,DELTA,x2_path,x3_path)

    x1_min = x1_(1);     x1_max = x1_(end);
    
    T = t_(end);             
    N = round(T/DELTA);      

    t_path = linspace(0,T,N+1);             

    x1_init = x0(1);
    x1_path = [x1_init zeros(1,N)];        
    
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

        [drift1,~,~] = b{i}(t_curr,x1_curr,x2_curr,x3_curr);
        
        dx1 = drift1 * DELTA; 

        x1_next = x1_curr + dx1;
        x1_next = min(max(x1_next, x1_min), x1_max);           % project x1 back onto the space grid 
        x1_path(n+1) = x1_next;
     
        prev_day = curr_day;

    end

end
