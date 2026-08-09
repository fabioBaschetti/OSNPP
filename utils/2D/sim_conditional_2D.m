function [t_path,x1_path,ttau_,iota_] = sim_conditional_2D(x0,i0,x1_,x2_,t_,p,b,DELTA,x2_path)

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

    for n = 1:N

        x1_curr = x1_path(n);
        x2_curr = x2_path(n);

        t_curr = t_path(n);
    
        [~,idx] = min(abs(t_(1:end-1) - t_curr));                         % index of the policy table to be used for time t_curr 

        action = 0;

        action = get_action_2D(p(:,:,:,idx),i,x1_curr,x2_curr,x1_,x2_);   % extract best action from current policy table
        if action ~= 0                                                    % if switch ...
            i = action;                                                   % ... update current regime
            ttau_ = horzcat(ttau_,t_curr);                                % ... update control
            iota_ = horzcat(iota_,i);                                     % "                "
        end

        [drift1, ~] = b{i}(t_curr,x1_curr,x2_curr);
            
        dx1 = drift1 * DELTA;

        x1_next = x1_curr + dx1;
        x1_next = min(max(x1_next, x1_min), x1_max);           % project x1 back onto the space grid 
        x1_path(n+1) = x1_next;

    end

end
