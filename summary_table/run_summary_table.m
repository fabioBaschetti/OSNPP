function results = run_summary_table()
%RUN_SUMMARY_TABLE Run existing scenario scripts and export summary metrics.
%
% The existing scripts are executed without modifying their code. Each script
% is expected to leave simulated paths such as t_path, x1_path, x2_path,
% ttau_, and iota_ in the MATLAB base workspace.

    this_dir = fileparts(mfilename('fullpath'));
    root_dir = fileparts(this_dir);
    out_dir  = fullfile(this_dir, 'output');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    scenarios = scenario_definitions(root_dir);
    results = table();
    original_dir = pwd;
    original_path = path;

    % Keep the post-processing helper functions visible even when the
    % existing scenario scripts change the current folder.
    addpath(this_dir);
    cleanupObj = onCleanup(@() restore_environment(original_dir, original_path)); %#ok<NASGU>

    for q = 1:numel(scenarios)
        sc = scenarios(q);
        fprintf('\n--- Running scenario %d/%d: %s ---\n', q, numel(scenarios), sc.Name);

        try
            evalin('base', 'clear variables');
            evalin('base', 'close all force');
            evalin('base', sprintf('cd(''%s'')', escape_single_quotes(sc.Folder)));
            evalin('base', sprintf('run(''%s'')', escape_single_quotes(fullfile(sc.Folder, sc.Script))));

            vars = collect_base_workspace();
            metrics = compute_summary_metrics(vars, sc.Name, sc.Economy);
            metrics.Status = "ok";
            metrics.ErrorMessage = "";
        catch ME
            warning('Scenario failed: %s\n%s', sc.Name, ME.message);
            metrics = empty_metric_row(sc.Name, sc.Economy);
            metrics.Status = "failed";
            metrics.ErrorMessage = string(ME.message);
        end

        results = [results; struct2table(metrics, 'AsArray', true)]; %#ok<AGROW>
    end

    write_summary_outputs(results, out_dir);
    fprintf('\nSummary outputs written to:\n%s\n', out_dir);
end

function s = escape_single_quotes(s)
    s = char(s);
    s = strrep(s, '''', '''''');
end

function row = empty_metric_row(name, economy)
    row = struct();
    row.Scenario = string(name);
    row.Economy = string(economy);
    row.HorizonDays = NaN;
    row.TotalCost = NaN;
    row.RunningCost = NaN;
    row.SwitchingCost = NaN;
    row.AvgAbsTrackingError = NaN;
    row.ShortageEnergy = NaN;
    row.ExcessEnergy = NaN;
    row.ShortageTimeSharePct = NaN;
    row.ExcessTimeSharePct = NaN;
    row.NumSwitches = NaN;
    row.TimeAtPminSharePct = NaN;
    row.TimeAtPmaxSharePct = NaN;
    row.PurchasesEnergy = NaN;
    row.SalesEnergy = NaN;
    row.PurchaseCost = NaN;
    row.SalesRevenue = NaN;
    row.AverageSellPrice = NaN;
    row.Status = "";
    row.ErrorMessage = "";
end

function restore_environment(original_dir, original_path)
    if exist(original_dir, 'dir')
        cd(original_dir);
    end
    path(original_path);
end
