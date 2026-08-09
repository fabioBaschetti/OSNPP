clear
close all
clc


cd("total_load\")

files = dir(fullfile(pwd, 'total_load_*'));
filenames = string({files.name}');
clear files

tab = table(datetime.empty(0,1),double.empty(0,1), ...
            'VariableNames', {'DateTime','Load'});

n = numel(filenames);
for i = 1:n

    filename = filenames(i);
    opts = detectImportOptions(filename);
    opts.SelectedVariableNames = opts.VariableNames(1:2);  
    opts.VariableNames(1:2) = {'DateTime', 'Load'};            
    tmp = readtable(filename, opts);
    tmp(end-1:end,:) = [];
    tmp = flipud(tmp);

    fn = char(filename);
    yy = str2double(fn(12:15));
    if double(yy) == 2025
        new_row = table(tmp.DateTime(1)-minutes(15),tmp.Load(1), 'VariableNames',tmp.Properties.VariableNames);
        tmp = [new_row; tmp];
    end
    clear fn yy new_row

    % ------ resolve start/end of daylight saving time
    idx = find(tmp.DateTime(2:end  ) - tmp.DateTime(1:end-1) > hours(1) + minutes(1));
    IDX = find(tmp.DateTime(2:end  ) - tmp.DateTime(1:end-1) < hours(0) + minutes(1));

    tmp.DateTime(idx+1:IDX(1)) = tmp.DateTime(idx+1:IDX(1)) - hours(1); 
    tmp.DateTime(IDX(2:4))     = tmp.DateTime(IDX(2:4))     - hours(1);
    tmp = sortrows(tmp, 1);
    % -----

    tab = [tab; tmp];

end
clear filename filenames i idx IDX n opts tmp

cd("..\")

isFeb29 = (month(tab.DateTime) == 2) & (day(tab.DateTime) == 29);
tab(isFeb29, :) = []; 
clear isFeb29

figure
plot(tab.DateTime,tab.Load)
grid on
ylabel('MW')
title('Total Load')

total_load = tab;
clear tab

save('total_load.mat','total_load')
