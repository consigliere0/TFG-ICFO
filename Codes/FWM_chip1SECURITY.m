% =========================================================================
%% TFG - COINCIDENCE COUNTS & CAR DATA PROCESSING (Standard Deviation)
% =========================================================================
clear; close all; format long g
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

tau = 1000 * 1e-12; 

% Paths and file listing
data_path = '../Tests/timetagger/wvg2/data0/';
files = dir(fullfile(data_path, '*2026-*'));

% Initialize arrays
powers_TE_mW = []; powers_TE_dBm = [];
powers_TM_mW = []; powers_TM_dBm = [];

rawCC_TE = []; rawCC_TE_std = []; acc_TE = []; acc_TE_std = [];
trueCC_TE = []; trueCC_TE_std = []; car_TE = []; %car_TE_sem = []; % CAR keeps SEM
car_TE_std = []; 

rawCC_TM = []; rawCC_TM_std = []; acc_TM = []; acc_TM_std = [];
trueCC_TM = []; trueCC_TM_std = []; car_TM = []; %car_TM_sem = [];
car_TM_std = [];

disp('Processing data files...');

for i = 1:length(files)
    fname = files(i).name;
    full_file_path = fullfile(files(i).folder, fname);
    
    token_TE = regexp(fname, 'TE(\d+)', 'tokens');
    token_TM = regexp(fname, 'TM(\d+)', 'tokens');
    
    if isempty(token_TE) && isempty(token_TM), continue; end
    
    data = readtable(full_file_path, 'FileType', 'text', 'Delimiter', ',', ...
        'NumHeaderLines', 9, 'VariableNamingRule', 'preserve');
    
    % Vectors
    CH1_vec = data{:, 3}; CH2_vec = data{:, 4}; Coinc_vec = data{:, 8};
    N = height(data);
    
    % Means & Std Devs
    CH1_m = mean(CH1_vec, 'omitnan'); CH1_s = std(CH1_vec, 'omitnan');
    CH2_m = mean(CH2_vec, 'omitnan'); CH2_s = std(CH2_vec, 'omitnan');
    raw_m = mean(Coinc_vec, 'omitnan'); raw_s = std(Coinc_vec, 'omitnan');
    
    % Accidentals
    acc_m = CH1_m * CH2_m * tau;
    acc_s = sqrt((CH2_m*tau*CH1_s)^2 + (CH1_m*tau*CH2_s)^2);
    
    % True Coincidences
    true_m = raw_m - acc_m;
    true_s = sqrt(raw_s^2 + acc_s^2);
    
    % CAR (Keep SEM here)
    car_m = true_m / acc_m;
    %car_sem = (car_m * sqrt((true_s/true_m)^2 + (acc_s/acc_m)^2)) / sqrt(N);
    car_std = (car_m * sqrt((true_s/true_m)^2 + (acc_s/acc_m)^2));
    
    if ~isempty(token_TE)
        p_dBm = str2double(token_TE{1}{1}) / 100;
        powers_TE_mW(end+1) = 10^(p_dBm/10); powers_TE_dBm(end+1) = p_dBm;
        rawCC_TE(end+1) = raw_m; rawCC_TE_std(end+1) = raw_s;
        acc_TE(end+1) = acc_m; acc_TE_std(end+1) = acc_s;
        trueCC_TE(end+1) = true_m; trueCC_TE_std(end+1) = true_s;
        car_TE(end+1) = car_m; car_TE_std(end+1) = car_std;
    else
        p_dBm = str2double(token_TM{1}{1}) / 100;
        powers_TM_mW(end+1) = 10^(p_dBm/10); powers_TM_dBm(end+1) = p_dBm;
        rawCC_TM(end+1) = raw_m; rawCC_TM_std(end+1) = raw_s;
        acc_TM(end+1) = acc_m; acc_TM_std(end+1) = acc_s;
        trueCC_TM(end+1) = true_m; trueCC_TM_std(end+1) = true_s;
        car_TM(end+1) = car_m; car_TM_std(end+1) = car_std;
    end
end

% Sort Arrays (Sorting logic omitted for brevity, use your previous method)
% [powers_TE_mW, idx] = sort(powers_TE_mW); ... etc

%% Plot 1: Raw CC and Accidentals (mW)
figure('Name', 'Raw CC', 'Color', 'w'); hold on; grid on; 
errorbar(powers_TE_mW, rawCC_TE, rawCC_TE_std, '-ob', 'LineWidth', 1, 'DisplayName', 'TE Raw CC');
errorbar(powers_TE_mW, acc_TE, acc_TE_std, '--sb', 'LineWidth', 1, 'DisplayName', 'TE Acc');
errorbar(powers_TM_mW, rawCC_TM, rawCC_TM_std, '-or', 'LineWidth', 1, 'DisplayName', 'TM Raw CC');
errorbar(powers_TM_mW, acc_TM, acc_TM_std, '--sr', 'LineWidth', 1, 'DisplayName', 'TM Acc');
xlabel('Laser Power (mW)'); ylabel('Counts (cps)'); legend;

%% Plot 2: True CC with Quadratic Fit (using Std Dev)
figure('Name', 'True CC', 'Color', 'w'); hold on; grid on;
errorbar(powers_TE_mW, trueCC_TE, trueCC_TE_std, 'ob', 'DisplayName', 'TE True CC');
errorbar(powers_TM_mW, trueCC_TM, trueCC_TM_std, 'or', 'DisplayName', 'TM True CC');
% Fit uses standard least squares (automatic weighting by variance if specified)
fit_te = fit(powers_TE_mW', trueCC_TE', 'poly2'); plot(fit_te, 'b');
fit_tm = fit(powers_TM_mW', trueCC_TM', 'poly2'); plot(fit_tm, 'r');
xlabel('Laser Power (mW)'); ylabel('True Coincidences (cps)');

%% Plot 3: CAR (using SEM)
figure('Name', 'CAR', 'Color', 'w'); hold on; grid on;
errorbar(powers_TE_mW, car_TE, car_TE_std, '-ob', 'DisplayName', 'TE CAR');
errorbar(powers_TM_mW, car_TM, car_TM_std, '-or', 'DisplayName', 'TM CAR');
xlabel('Laser Power (mW)'); ylabel('CAR'); legend;