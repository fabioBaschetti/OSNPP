function row = append_current_workspace_to_summary(scenario_name, economy)
%APPEND_CURRENT_WORKSPACE_TO_SUMMARY Summarize the current base workspace.
%
% Usage:
%   row = append_current_workspace_to_summary('Closed benchmark','closed')
%   row = append_current_workspace_to_summary('Open economy, price maker','open')
%
% This is useful when you prefer to run the existing scenario scripts manually.

    if nargin < 2
        error('Provide both scenario_name and economy, e.g. ''closed'' or ''open''.');
    end

    this_dir = fileparts(mfilename('fullpath'));
    out_dir  = fullfile(this_dir, 'output');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    vars = collect_base_workspace();
    row_struct = compute_summary_metrics(vars, string(scenario_name), string(economy));
    row_struct.Status = "ok";
    row_struct.ErrorMessage = "";
    row = struct2table(row_struct, 'AsArray', true);

    csv_file = fullfile(out_dir, 'summary_table_full.csv');
    if exist(csv_file, 'file')
        old = readtable(csv_file, 'TextType', 'string');
        results = [old; row];
    else
        results = row;
    end

    write_summary_outputs(results, out_dir);
end
