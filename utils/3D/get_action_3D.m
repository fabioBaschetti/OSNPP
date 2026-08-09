function out = get_action_3D(p_slice,i,y1,y2,y3,x1_,x2_,x3_)
    k1 = get_nearest_grid_node(y1,x1_);
    k2 = get_nearest_grid_node(y2,x2_);
    k3 = get_nearest_grid_node(y3,x3_);
    out = p_slice(k2,k1,k3,i);
end

function out = get_nearest_grid_node(y,x_)
    Nx = length(x_);
    x0 = x_(1);
    dx = (x_(end) - x0) / (Nx - 1);
    out = max(1, min(Nx, round(1+(y-x0)/dx)));
end
