clear
close all
clc


% ----- get baseline capacity for each RES as of 31/12/2023
cd("renewables\")

filename = 'installed_capacity_at_31122023.xlsx';
opts = detectImportOptions(filename);
opts.SelectedVariableNames = opts.VariableNames([7,9]);  
opts.VariableNames([7,9]) = {'Source', 'MW'};            
tab = readtable(filename, opts);
tab(end-1:end,:) = [];

clear filename opts

cd("..\")

groups = groupsummary(tab, 'Source', 'sum', 'MW');
clear tab
groups = groups(:, {'Source', 'sum_MW'});
groups.Properties.VariableNames = {'Source', 'TOT_Capacity'};

baseline_capacity = struct();

for i = 1:height(groups)
    src = groups.Source{i};
    if ~isempty(src)
        name = matlab.lang.makeValidName(src);
        baseline_capacity.(name) = groups.TOT_Capacity(i);
    end
end
clear i src name groups
% -----


% ----- get running  capacity for each RES between 01/01/2024 and 31/12/2025
%       (with monthly resolution)
cd("renewables\")

files = dir(fullfile(pwd, 'installed_capacity_in_year*'));                     %
filenames = string({files.name}');
clear files
filename = filenames(1);
opts = detectImportOptions(filename);
opts.SelectedVariableNames = opts.VariableNames([1,2,7,9]);  
opts.VariableNames([1,2,7,9]) = {'Year', 'Month', 'Source', 'MW'};  
tab = readtable(filename, opts);                                               % new installed capacity over the years    
tab(end-1:end,:) = [];                                                         %
for i = 2:numel(filenames)                                                     %
    filename = filenames(i);                                                   %
    tmp = readtable(filename, opts);                                           %
    tmp(end-1:end,:) = [];                                                     %
    tab = [tab; tmp];                                                          %
end
clear i filename tmp
tab.Year = str2double(tab.Year);

clear filenames opts

cd("..\")

groups = groupsummary(tab, {'Source', 'Year', 'Month'}, 'sum', 'MW');          % group (i.e. sum) new installations based on (Source,Year,Month)
clear tab
groups = groups(:, {'Source', 'Year','Month', 'sum_MW'});

yy = table2array(unique(groups(:, 'Year')));
mm = (1:12)';
YY = repelem(yy, numel(mm),1);
MM = repmat(mm, numel(yy),1);
YM = table(YY,MM, 'VariableNames',{'Year','Month'});                           % (°)
clear yy YY mm MM

RES = table(fieldnames(baseline_capacity), 'VariableNames',{'Source'});        % (*)

groups = outerjoin([repelem(RES,height(YM),1), repmat(YM,height(RES),1)], ...  % fill in the gaps: every RES (*) and YM (°) pair gets its own piece
                   groups, ...                                                 %                   of installed capacity ...
                   'Keys', {'Source','Year','Month'}, ...
                   'MergeKeys', true, ...
                   'Type', 'left');
if any(ismissing(groups.sum_MW)) 
    groups.sum_MW(ismissing(groups.sum_MW)) = 0;                               % ... eventually zero (if no installation  on (Source,Year,Month))  
end
clear YM RES

groups.cumsum_MW = grouptransform(groups, "Source", ...                        % add column 'cumsum_MW' for the comulative sum of installed capacity
                                          @(x) cumsum(x), "sum_MW").sum_MW;    % over [01/01/2024,t] (for all RES)

groups.installed_capacity_run = zeros(height(groups),1);                       % add baseline capacity at 31/12/2023
for i = 1:height(groups)
    src = groups.Source{i};
    groups.installed_capacity_run(i) = groups.cumsum_MW(i) + baseline_capacity.(matlab.lang.makeValidName(src));
end
clear i src

groups = groups(:, {'Source', 'Year','Month', 'installed_capacity_run'});

clear baseline_capacity

max_cap = groupsummary(groups, 'Source', 'max', 'installed_capacity_run');     % get end-of-period (eop) capacity 
max_cap = max_cap(:, {'Source', 'max_installed_capacity_run'});
max_cap.Properties.VariableNames = {'Source', 'installed_capacity_eop'};

groups = outerjoin(groups, max_cap(:, {'Source','installed_capacity_eop'}),... % join eop capacity as a new column
                   'Keys', 'Source', ...
                   'MergeKeys', true, ...
                   'Type', 'left');

clear max_cap

mapNames = containers.Map({'Bioenergie','Geotermoelettrico','Idroelettrico','Solare'      ,'Eolico'}, ...
                          {'Biomass'   ,'Geothermal'       ,'Hydro'        ,'Photovoltaic','Wind'  });
for i = 1:height(groups)
    if isKey(mapNames, groups.Source{i})
        groups.Source{i} = mapNames(groups.Source{i});
    end
end
clear i mapNames

Y = groups.Year ;                                                              % replace (Y,M) with Start-End datetime objects
M = groups.Month;
Start = datetime(Y,M,1,0,0,0);
End   = dateshift(Start, 'end', 'month') + duration(23,59,59);
Start.Format = 'dd-MMM-yyyy HH:mm:ss';
clear Y M

tot_capacity_RES = table( ...
    groups.Source, ...
    Start, ...
    End, ...
    groups.installed_capacity_run, ...
    groups.installed_capacity_eop, ...
    'VariableNames', {'Source','Start','End','installed_capacity_run','installed_capacity_eop'} );

clear groups Start End
% -----


% ----- time serie of renewable production between 01/01/2024 and 31/12/2025
cd("renewables\")

files = dir(fullfile(pwd, 'production_*'));
filenames = string({files.name}');
clear files

n = numel(filenames);
for i = 1:n

    filename = filenames(i);
    opts = detectImportOptions(filename);
    opts.SelectedVariableNames = opts.VariableNames([1,3,2]);
    opts.VariableNames([1,3,2]) = {'DateTime','Source','Production'};
    tab = readtable(filename, opts);
    tab(end-1:end,:) = [];
    tab = flipud(tab);

    fn = char(filename);
    yy = str2double(fn(12:15));
    if yy<=2024
        DELTA = hours(1);              % renewable production was measured at the hourly scale in 2024
    else
        DELTA = minutes(15);           % but they eventually switched to 15-minutes   starting in 2025
    end
    clear fn yy

    sources = unique(tab.Source);

    tmpsBySource = struct();
    for k = 1:numel(sources)
        
        src = sources(k);
        src = src{1};
        rows = string(tab.Source) == src;
        tmpsBySource.(matlab.lang.makeValidName(char(src))) = tab(rows,[1,3]);

        if strcmp(src,'Biomass') == 1       % production at 00.00 on January 1 is missing for Biomass 
            new_row = table(tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(1)-DELTA,tmpsBySource.(matlab.lang.makeValidName(char(src))).Production(1), ...
                            'VariableNames',tmpsBySource.(matlab.lang.makeValidName(char(src))).Properties.VariableNames);
            tmpsBySource.(matlab.lang.makeValidName(char(src))) = [new_row; tmpsBySource.(matlab.lang.makeValidName(char(src)))]; 
            clear new_row
        end

        % ------ resolve start/end of daylight saving time
        idx = find(tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(2:end  ) - ...
                   tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(1:end-1) > ...
                   hours(1) + minutes(1));

        IDX = find(tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(2:end  ) - ...
                   tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(1:end-1) < ...
                   hours(0) + minutes(1));

        tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(idx+1:IDX) = tmpsBySource.(matlab.lang.makeValidName(char(src))).DateTime(idx+1:IDX) - hours(1); 
        % ------

    end

    if i == 1
        tabsBySource = tmpsBySource;
    else
        for k = 1:numel(sources)
            src = sources(k);
            src = src{1};
            tabsBySource.(matlab.lang.makeValidName(char(src))) = [tabsBySource.(matlab.lang.makeValidName(char(src))); ...
                                                                   tmpsBySource.(matlab.lang.makeValidName(char(src)))];
        end
    end

end
clear DELTA filename filenames i idx IDX k n opts rows src tab tmpsBySource

cd("..\")


% ----- rescale renewable production at eop capacity
for i = 1:numel(sources)               
                                                                               
    source = sources{i};
    tmp = tot_capacity_RES(string(tot_capacity_RES.Source)==matlab.lang.makeValidName(source),:);

    ddtt = tabsBySource.(matlab.lang.makeValidName(source)).DateTime;
    for n = 1:numel(ddtt)
        idx = find((ddtt(n)>=tmp.Start) & (ddtt(n)<=tmp.End));
        tabsBySource.(matlab.lang.makeValidName(source)).Production(n) = 1000 ...
                                                                         *(tmp.installed_capacity_eop(idx)/tmp.installed_capacity_run(idx)) ...
                                                                         * tabsBySource.(matlab.lang.makeValidName(source)).Production(n);   
    end

end
clear i source tmp ddtt n idx tot_capacity_RES
% -----


% ----- interpolate production date to finest scale (i.e. 15 minutes)
for k = 1:numel(sources)

    src = sources(k);
    src = src{1};
    tt = tabsBySource.(matlab.lang.makeValidName(char(src)));

    ts = timetable(tt.DateTime,tt.Production, 'VariableNames',{'Production'});

    ts = retime(ts, 'regular', 'linear', 'TimeStep', minutes(15));             

    tt = timetable2table(ts, 'ConvertRowTimes',true);
    tt.Properties.VariableNames{'Time'} = 'DateTime';

    tabsBySource.(matlab.lang.makeValidName(char(src))) = tt;

end
clear k src ts tt
% -----

% ----- remove Feb 29
for k = 1:numel(sources)
    
    src = sources(k);
    src = src{1};

    isFeb29 = (month(tabsBySource.(matlab.lang.makeValidName(char(src))).DateTime) == 2) & (day(tabsBySource.(matlab.lang.makeValidName(char(src))).DateTime) == 29); 

    tabsBySource.(matlab.lang.makeValidName(char(src)))(isFeb29, :) = [];

end
clear k src isFeb29 
% -----


for k = 1:numel(sources)

    src = sources(k);
    src = src{1};

    figure
    plot(tabsBySource.(matlab.lang.makeValidName(char(src))).DateTime, tabsBySource.(matlab.lang.makeValidName(char(src))).Production)
    grid on
    ylabel('MW')
    title(matlab.lang.makeValidName(char(src)))
    
end
clear k src

clear sources


fields = fieldnames(tabsBySource);
fields = fields(3:5);        % remove biomass and geothermic because of wierd graphs

X = cellfun(@(f) tabsBySource.(f).Production, fields, 'UniformOutput', false);
X = cat(2, X{:});

tot = sum(X, 2);
clear X

renewables = table(tabsBySource.(matlab.lang.makeValidName(fields{1})).DateTime, tot, 'VariableNames',{'DateTime','Production'});
clear fields tot

figure
plot(renewables.DateTime, renewables.Production)
grid on
ylabel('MW')
title('renewables')

save('renewables.mat','renewables')
