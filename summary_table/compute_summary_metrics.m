function row = compute_summary_metrics(vars, scenario_name, economy)
%COMPUTE_SUMMARY_METRICS Compute economic metrics from simulated paths.
%
% Required fields in vars:
%   t_path, x1_path, x2_path
% Optional fields:
%   x3_path, z_path, S_path, ttau_, iota_, i0, G, g, lambda1, lambda2,
%   gamma1, delta, P_min/P_max or P_min_IT/P_max_IT.

    require_field(vars, 't_path');
    require_field(vars, 'x1_path');
    require_field(vars, 'x2_path');

    t = rowvec(vars.t_path);
    P = rowvec(vars.x1_path);
    Y = rowvec(vars.x2_path);

    if numel(t) ~= numel(P) || numel(t) ~= numel(Y)
        error('t_path, x1_path, and x2_path must have the same length.');
    end
    if numel(t) < 2
        error('The simulated path must contain at least two time points.');
    end

    dt_vec = diff(t);
    T_horizon = t(end) - t(1);

    P_left = P(1:end-1);
    Y_left = Y(1:end-1);
    shortage = max(Y_left - P_left, 0);
    excess   = max(P_left - Y_left, 0);

    shortage_energy = sum(dt_vec .* shortage);
    excess_energy   = sum(dt_vec .* excess);

    avg_abs_tracking = sum(dt_vec .* abs(P_left - Y_left)) / T_horizon;
    shortage_share   = 100 * sum(dt_vec .* (Y_left > P_left)) / T_horizon;
    excess_share     = 100 * sum(dt_vec .* (P_left > Y_left)) / T_horizon;

    [P_min, P_max] = get_production_bounds(vars, P);
    tol = get_bound_tolerance(vars, P_min, P_max);
    time_at_pmin = 100 * sum(dt_vec .* (abs(P_left - P_min) <= tol)) / T_horizon;
    time_at_pmax = 100 * sum(dt_vec .* (abs(P_left - P_max) <= tol)) / T_horizon;

    [running_cost, purchase_cost, sales_revenue, avg_sell_price] = ...
        compute_running_cost(vars, economy, dt_vec, P_left, Y_left, shortage, excess);

    switching_cost = compute_switching_cost(vars);
    total_cost = running_cost + switching_cost;

    is_open = strcmpi(string(economy), "open");
    if is_open
        purchases_energy = shortage_energy;
        sales_energy = excess_energy;
    else
        purchases_energy = NaN;
        sales_energy = NaN;
        purchase_cost = NaN;
        sales_revenue = NaN;
        avg_sell_price = NaN;
    end

    row = struct();
    row.Scenario = string(scenario_name);
    row.Economy = string(economy);
    row.HorizonDays = T_horizon;
    row.TotalCost = total_cost;
    row.RunningCost = running_cost;
    row.SwitchingCost = switching_cost;
    row.AvgAbsTrackingError = avg_abs_tracking;
    row.ShortageEnergy = shortage_energy;
    row.ExcessEnergy = excess_energy;
    row.ShortageTimeSharePct = shortage_share;
    row.ExcessTimeSharePct = excess_share;
    row.NumSwitches = get_num_switches(vars);
    row.TimeAtPminSharePct = time_at_pmin;
    row.TimeAtPmaxSharePct = time_at_pmax;
    row.PurchasesEnergy = purchases_energy;
    row.SalesEnergy = sales_energy;
    row.PurchaseCost = purchase_cost;
    row.SalesRevenue = sales_revenue;
    row.AverageSellPrice = avg_sell_price;
end

function require_field(vars, nm)
    if ~isfield(vars, nm)
        error('Missing required variable: %s.', nm);
    end
end

function x = rowvec(x)
    x = x(:).';
end

function [P_min, P_max] = get_production_bounds(vars, P)
    if isfield(vars, 'P_min_IT')
        P_min = vars.P_min_IT;
    elseif isfield(vars, 'P_min')
        P_min = vars.P_min;
    else
        P_min = min(P);
    end

    if isfield(vars, 'P_max_IT')
        P_max = vars.P_max_IT;
    elseif isfield(vars, 'P_max')
        P_max = vars.P_max;
    else
        P_max = max(P);
    end
end

function tol = get_bound_tolerance(vars, P_min, P_max)
    if isfield(vars, 'x1_') && numel(vars.x1_) > 1
        xgrid = rowvec(vars.x1_);
        tol = 0.5 * min(diff(xgrid));
    else
        tol = max(1e-8, 1e-4 * max(1, P_max - P_min));
    end
end

function [running_cost, purchase_cost, sales_revenue, avg_sell_price] = compute_running_cost(vars, economy, dt_vec, P, Y, shortage, excess)
    gamma1 = get_scalar(vars, 'gamma1', 0);
    is_open = strcmpi(string(economy), "open");

    purchase_cost = NaN;
    sales_revenue = NaN;
    avg_sell_price = NaN;

    if is_open
        S = get_sell_price_path(vars, numel(P));
        delta = get_scalar(vars, 'delta', 0);
        purchase_cost = sum(dt_vec .* (S + delta) .* shortage);
        sales_revenue = sum(dt_vec .* S .* excess);
        running_cost = purchase_cost - sales_revenue + sum(dt_vec .* gamma1 .* P);
        avg_sell_price = sum(dt_vec .* S) / sum(dt_vec);
    else
        lambda1 = get_scalar(vars, 'lambda1', NaN);
        lambda2 = get_scalar(vars, 'lambda2', NaN);
        if isnan(lambda1) || isnan(lambda2)
            running_cost = compute_running_cost_from_fcell(vars, dt_vec, P, Y);
        else
            running_cost = sum(dt_vec .* (lambda1 .* excess + lambda2 .* shortage + gamma1 .* P));
        end
    end
end

function S = get_sell_price_path(vars, n_left)
    if isfield(vars, 'S_path')
        S_full = rowvec(vars.S_path);
        S = S_full(1:n_left);
        return;
    end

    if isfield(vars, 'x3_path') && isfield(vars, 'P_max_EU') && ...
            isfield(vars, 's_res') && isfield(vars, 's_npp') && isfield(vars, 's_fss')
        x3 = rowvec(vars.x3_path);
        x3 = x3(1:n_left);
        S = vars.s_npp * ones(size(x3));
        S(x3 <= 0) = vars.s_res;
        S(x3 > vars.P_max_EU) = vars.s_fss;
        return;
    end

    error('Open-economy metrics require S_path or enough variables to reconstruct it.');
end

function out = compute_running_cost_from_fcell(vars, dt_vec, P, Y)
    if ~isfield(vars, 'f') || ~iscell(vars.f)
        error('Cannot compute running cost: missing lambda parameters and f cell array.');
    end
    regimes = reconstruct_regime_path(vars, numel(P));
    out = 0;
    for n = 1:numel(P)
        i = regimes(n);
        out = out + dt_vec(n) * vars.f{i}(0, P(n), Y(n));
    end
end

function regimes = reconstruct_regime_path(vars, n_left)
    i0 = get_scalar(vars, 'i0', 1);
    regimes = i0 * ones(1, n_left);
    if isfield(vars, 'ttau_') && isfield(vars, 'iota_') && isfield(vars, 't_path')
        t = rowvec(vars.t_path);
        ttau = rowvec(vars.ttau_);
        iota = rowvec(vars.iota_);
        for k = 1:min(numel(ttau), numel(iota))
            regimes(t(1:n_left) >= ttau(k)) = iota(k);
        end
    end
end

function switching_cost = compute_switching_cost(vars)
    switching_cost = 0;
    if ~isfield(vars, 'ttau_') || ~isfield(vars, 'iota_') || isempty(vars.ttau_)
        return;
    end

    ttau = rowvec(vars.ttau_);
    iota = rowvec(vars.iota_);
    i0 = get_scalar(vars, 'i0', NaN);
    if isnan(i0)
        error('Switches are present but i0 is missing.');
    end

    prev = [i0, iota(1:end-1)];
    next = iota;

    for k = 1:numel(next)
        i = prev(k);
        j = next(k);
        t_switch = ttau(min(k, numel(ttau)));
        switching_cost = switching_cost + single_switch_cost(vars, i, j, t_switch);
    end
end

function c = single_switch_cost(vars, i, j, t_switch)
    if isfield(vars, 'g') && isa(vars.g, 'function_handle')
        try
            c = vars.g(i, j, t_switch, 0, 0);
            return;
        catch
            % Fall back to G below.
        end
    end

    if isfield(vars, 'G')
        c = vars.G(i, j);
        return;
    end

    c = 0;
end

function n = get_num_switches(vars)
    if isfield(vars, 'ttau_')
        n = numel(vars.ttau_);
    else
        n = 0;
    end
end

function val = get_scalar(vars, nm, default_val)
    if isfield(vars, nm) && isscalar(vars.(nm))
        val = vars.(nm);
    else
        val = default_val;
    end
end
