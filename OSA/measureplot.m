clc; close all; clear all;

%{
data1 = 'W0000.CSV';

options = detectImportOptions(data1);
options.DataLine = [30, Inf];
options.VariableNamingRule = 'preserve';

datatable = readtable(data1, options);
disp('Showing first entries');

lam = datatable{:, 1};
pow = datatable{:,2};

figure;
plot(lam, pow, '-k')
%}

baseName = 'W000';
options = detectImportOptions([baseName, '0.CSV']);
options.DataLine = [30, Inf];
options.VariableNamingRule = 'preserve';

figure; hold on; grid on; box on;

numFile = 0:7;

for i = numFile
    nameFile = sprintf('%s%d.CSV', baseName, i);
    try
        data = readtable(nameFile, options);
        wvl = data{:, 1};
        pwr = data{:, 2};
        plot(wvl, pwr, sprintf('File W000%d', i))

    catch
        fprintf("Couldn't load file %s\n", nameFile);
    end
end