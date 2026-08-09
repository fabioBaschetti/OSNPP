function m = my_interp3D(M,grid,x1_out,x2_out,x3_out)

    % linear interpolation over a 3-dim lattice \mathcal{X} of size (Nx2 x Nx1 x Nx3)
    % extrapolation: clamping
    % 
    % (achtung! the order is the same as the one implied by matlab meshgrid built-in function)

    Nx1 = grid.Nx1;
    Nx2 = grid.Nx2;
    Nx3 = grid.Nx3;

    % get fractional indices
    tau1 = 1 + (x1_out - grid.x1_offset) * grid.inv_h1;
    tau2 = 1 + (x2_out - grid.x2_offset) * grid.inv_h2;
    tau3 = 1 + (x3_out - grid.x3_offset) * grid.inv_h3;

    % CLAMP outside borders
    tau1(tau1<=  1) =   1 + 1e-12;
    tau1(tau1>=Nx1) = Nx1 - 1e-12;
    tau2(tau2<=  1) =   1 + 1e-12;
    tau2(tau2>=Nx2) = Nx2 - 1e-12;
    tau3(tau3<=  1) =   1 + 1e-12;
    tau3(tau3>=Nx3) = Nx3 - 1e-12;

    i = floor(tau1);    % always in [1,Nx1-1]
    j = floor(tau2);    % always in [1,Nx2-1]
    k = floor(tau3);    % always in [1,Nx3-1]

    % weights
    w1 = tau1 - i;      % always in [0,1]
    w2 = tau2 - j;      % always in [0,1]
    w3 = tau3 - k;      % always in [0,1]

    u1 = 1 - w1;
    u2 = 1 - w2;
    u3 = 1 - w3;

    % implement explicit linear index formula
    % idx(j,i,k) = j     + (i-1)*Nx2 + (k-1)*Nx2*Nx1
    % == base == 
    %                s_j         s_i         - s_k - 

    shift_i = Nx2;
    shift_j = 1;
    shift_k = Nx2 * Nx1;

    base = j + (i-1)*Nx2 + (k-1)*shift_k;

    M000 = M(base);
    M100 = M(base + shift_i);
    M010 = M(base + shift_j);
    M110 = M(base + shift_i + shift_j);
    M001 = M(base + shift_k);
    M101 = M(base + shift_i + shift_k);
    M011 = M(base + shift_j + shift_k);
    M111 = M(base + shift_i + shift_j + shift_k);

    m = u1.*u2.*u3.*M000 + ...
        w1.*u2.*u3.*M100 + ...
        u1.*w2.*u3.*M010 + ...
        w1.*w2.*u3.*M110 + ...
        u1.*u2.*w3.*M001 + ...
        w1.*u2.*w3.*M101 + ...
        u1.*w2.*w3.*M011 + ...
        w1.*w2.*w3.*M111;

end
