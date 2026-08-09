function scenarios = scenario_definitions(root_dir)
%SCENARIO_DEFINITIONS List of existing scripts to run for the summary table.
%
% This file only points to existing scripts. It does not change the numerical
% algorithm or the implementation of the solver.

    scenarios = struct([]);

    scenarios(end+1).Name    = "Closed benchmark";
    scenarios(end).Economy   = "closed";
    scenarios(end).Folder    = fullfile(root_dir, 'OSNPP_clsdmkt');
    scenarios(end).Script    = 'main.m';

    scenarios(end+1).Name    = "Open economy, price taker";
    scenarios(end).Economy   = "open";
    scenarios(end).Folder    = fullfile(root_dir, 'OSNPP_openmkt_pricetaker');
    scenarios(end).Script    = 'main.m';

    scenarios(end+1).Name    = "Open economy, price maker";
    scenarios(end).Economy   = "open";
    scenarios(end).Folder    = fullfile(root_dir, 'OSNPP_openmkt_pricemaker');
    scenarios(end).Script    = 'main.m';
end
