function m = my_interp2D(M,x1_in,x2_in,x1_out,x2_out,extrap)

    % linear interpolation over a 2-dim lattice  \mathcal{X} \in (Nx2 x Nx1) 
    % + extrapolation: clamping
    %                  linear
    % in particular \mathcal{X} \equiv [X1_in,X2_in] = meshgrid(x1_in,x2_in)

    Nx1 = length(x1_in);
    Nx2 = length(x2_in);

    h1 = (x1_in(end) - x1_in(1)) / (Nx1 - 1);    % (uniform) step size on the x1_in grid
    h2 = (x2_in(end) - x2_in(1)) / (Nx2 - 1);    % (uniform) step size on the x2_in grid

    % get fractional indices
    tau1 = 1 + (x1_out - x1_in(1)) / h1;
    tau2 = 1 + (x2_out - x2_in(1)) / h2;

    switch lower(extrap)

        case 'clamp'
            % and clamp outside borders
            tau1(tau1<=  1) =   1 + 1e-12;
            tau1(tau1>=Nx1) = Nx1 - 1e-12;
            tau2(tau2<=  1) =   1 + 1e-12;
            tau2(tau2>=Nx2) = Nx2 - 1e-12;

            i = floor(tau1);    % always in [1,Nx1-1]
            j = floor(tau2);    % always in [1,Nx2-1]

            w1 = tau1 - i;
            w2 = tau2 - j;

        case 'linear'

            i = min(max(floor(tau1), 1), Nx1-1);
            j = min(max(floor(tau2), 1), Nx2-1);

            w1 = tau1 - i;      % can be <0 or >1
            w2 = tau2 - j;      % can be <0 or >1

        otherwise
            error('ERROR: invalid mode');

    end

    idx00 = sub2ind([Nx2,Nx1],j  ,i  );
    idx10 = sub2ind([Nx2,Nx1],j  ,i+1);
    idx01 = sub2ind([Nx2,Nx1],j+1,i  );
    idx11 = sub2ind([Nx2,Nx1],j+1,i+1);

    M00 = M(idx00);
    M10 = M(idx10);
    M01 = M(idx01);
    M11 = M(idx11);

    m = (1-w1).*(1-w2).*M00 + ...
           w1 .*(1-w2).*M10 + ...
        (1-w1).*   w2 .*M01 + ...
           w1 .*   w2 .*M11;

end
