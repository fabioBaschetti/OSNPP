function out = get_action_2D(p,i,y1,y2,x1_,x2_)
    k1 = get_nearest_grid_node(y1,x1_);
    k2 = get_nearest_grid_node(y2,x2_);
    out = p(k2,k1,i);
end

function out = get_nearest_grid_node(y,x_)
    Nx = length(x_);
    x0 = x_(1);
    dx = (x_(end) - x0) / (Nx - 1);
    out = max(1, min(Nx, round(1+(y-x0)/dx)));
end
