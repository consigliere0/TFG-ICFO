% =========================================================================
%% TFG - COINCIDENCE COUNTS & CAR DATA PROCESSING
% =========================================================================
%%                     SECTION 1: INITIAL SETTINGS
% =========================================================================
clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Coincidence window = 1000 ps (from file header)
tau = 1000 * 1e-12;

% =========================================================================
%% SECTION 2: LOAD DATA
% =========================================================================
data_path = '../Tests/timetagger/wvg2/data1/';
files = dir(fullfile(data_path, '*2026-*'));

% Initialize arrays (mW and dBm stored in parallel)
powers_TE_mW  = []; powers_TE_dBm = [];
powers_TM_mW  = []; powers_TM_dBm = [];
rawCC_TE = []; rawCC_TE_std = []; acc_TE = []; acc_TE_std = [];
trueCC_TE = []; trueCC_TE_std = []; car_TE = []; car_TE_std = [];
rawCC_TM = []; rawCC_TM_std = []; acc_TM = []; acc_TM_std = [];
trueCC_TM = []; trueCC_TM_std = []; car_TM = []; car_TM_std = [];

disp('Processing data files...');

for i = 1:length(files)
    fname = files(i).name;
    full_file_path = fullfile(files(i).folder, fname);

    token_TE = regexp(fname, 'TE(\d+)', 'tokens');
    token_TM = regexp(fname, 'TM(\d+)', 'tokens');

    if isempty(token_TE) && isempty(token_TM)
        fprintf('Skipping (no TE/TM token found): %s\n', fname);
        continue;
    end

    try
        data = readtable(full_file_path, 'FileType', 'text', 'Delimiter', ',', ...
            'NumHeaderLines', 9, 'VariableNamingRule', 'preserve');
    catch
        warning('Could not read file: %s', full_file_path);
        continue;
    end

    if size(data, 2) < 8
        warning('File %s has fewer than 8 columns after parsing. Skipping...', fname);
        continue;
    end

    % =========================================================================
    %% SECTION 3: MAGNITUDE CALCULATION
    % =========================================================================
    CH1_vec   = data{:, 3};
    CH2_vec   = data{:, 4};
    Coinc_vec = data{:, 8};

    % Means
    CH1_mean       = mean(CH1_vec,   'omitnan');
    CH2_mean       = mean(CH2_vec,   'omitnan');
    Coinc_raw_mean = mean(Coinc_vec, 'omitnan');

    % Standard deviations (sample std, N-1)
    CH1_std       = std(CH1_vec,   0, 'omitnan');
    CH2_std       = std(CH2_vec,   0, 'omitnan');
    Coinc_raw_std = std(Coinc_vec, 0, 'omitnan');

    % Accidentals: A = CH1 * CH2 * tau
    accidentals     = CH1_mean * CH2_mean * tau;
    % Propagation: sigma_A = tau * sqrt((CH2*sigma_CH1)^2 + (CH1*sigma_CH2)^2)
    accidentals_std = sqrt((CH2_mean * tau * CH1_std)^2 + (CH1_mean * tau * CH2_std)^2);

    % True coincidences: C_true = C_raw - A
    true_coincidences     = Coinc_raw_mean - accidentals;
    % Propagation: sigma_true = sqrt(sigma_raw^2 + sigma_A^2)
    true_coincidences_std = sqrt(Coinc_raw_std^2 + accidentals_std^2);

    % CAR: C_true / A
    CAR     = true_coincidences / accidentals;
    % Propagation: sigma_CAR = CAR * sqrt((sigma_true/C_true)^2 + (sigma_A/A)^2)
    CAR_std = abs(CAR) * sqrt((true_coincidences_std / true_coincidences)^2 + (accidentals_std / accidentals)^2);

    if ~isempty(token_TE)
        power_dBm = str2double(token_TE{1}{1}) / 100;
        power_mW  = 10^(power_dBm / 10);
        powers_TE_mW(end+1)   = power_mW;
        powers_TE_dBm(end+1)  = power_dBm;
        rawCC_TE(end+1)       = Coinc_raw_mean; rawCC_TE_std(end+1)   = Coinc_raw_std;
        acc_TE(end+1)         = accidentals;    acc_TE_std(end+1)     = accidentals_std;
        trueCC_TE(end+1)      = true_coincidences; trueCC_TE_std(end+1) = true_coincidences_std;
        car_TE(end+1)         = CAR;            car_TE_std(end+1)     = CAR_std;
        fprintf('TE file: %s | %.2f dBm = %.2f mW | RawCC=%.2f | Acc=%.4f | TrueCC=%.2f | CAR=%.2f\n', ...
            fname, power_dBm, power_mW, Coinc_raw_mean, accidentals, true_coincidences, CAR);
    elseif ~isempty(token_TM)
        power_dBm = str2double(token_TM{1}{1}) / 100;
        power_mW  = 10^(power_dBm / 10);
        powers_TM_mW(end+1)   = power_mW;
        powers_TM_dBm(end+1)  = power_dBm;
        rawCC_TM(end+1)       = Coinc_raw_mean; rawCC_TM_std(end+1)   = Coinc_raw_std;
        acc_TM(end+1)         = accidentals;    acc_TM_std(end+1)     = accidentals_std;
        trueCC_TM(end+1)      = true_coincidences; trueCC_TM_std(end+1) = true_coincidences_std;
        car_TM(end+1)         = CAR;            car_TM_std(end+1)     = CAR_std;
        fprintf('TM file: %s | %.2f dBm = %.2f mW | RawCC=%.2f | Acc=%.4f | TrueCC=%.2f | CAR=%.2f\n', ...
            fname, power_dBm, power_mW, Coinc_raw_mean, accidentals, true_coincidences, CAR);
    end
end

fprintf('\nFiles parsed: %d TE, %d TM\n', length(powers_TE_mW), length(powers_TM_mW));

% =========================================================================
%% SECTION 4: SORT
% =========================================================================
[powers_TE_mW, sortIdx_TE] = sort(powers_TE_mW);
powers_TE_dBm    = powers_TE_dBm(sortIdx_TE);
rawCC_TE         = rawCC_TE(sortIdx_TE);   rawCC_TE_std   = rawCC_TE_std(sortIdx_TE);
acc_TE           = acc_TE(sortIdx_TE);     acc_TE_std     = acc_TE_std(sortIdx_TE);
trueCC_TE        = trueCC_TE(sortIdx_TE);  trueCC_TE_std  = trueCC_TE_std(sortIdx_TE);
car_TE           = car_TE(sortIdx_TE);     car_TE_std     = car_TE_std(sortIdx_TE);

[powers_TM_mW, sortIdx_TM] = sort(powers_TM_mW);
powers_TM_dBm    = powers_TM_dBm(sortIdx_TM);
rawCC_TM         = rawCC_TM(sortIdx_TM);   rawCC_TM_std   = rawCC_TM_std(sortIdx_TM);
acc_TM           = acc_TM(sortIdx_TM);     acc_TM_std     = acc_TM_std(sortIdx_TM);
trueCC_TM        = trueCC_TM(sortIdx_TM);  trueCC_TM_std  = trueCC_TM_std(sortIdx_TM);
car_TM           = car_TM(sortIdx_TM);     car_TM_std     = car_TM_std(sortIdx_TM);

disp('Data successfully processed. Generating plots...');

% =========================================================================
%% SECTION 5: PLOTS — LINEAR POWER AXIS (mW)
% =========================================================================

%% --- PLOT 1: Raw Coincidences & Accidentals vs Power (mW) ---
figure('Name', 'Raw CC and Accidentals vs Power (mW)', 'Color', 'w');
hold on; box on; grid on;

errorbar(powers_TE_mW, rawCC_TE, rawCC_TE_std, '-ob', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'CapSize', 6, 'DisplayName', 'TE Raw Coincidences');
errorbar(powers_TE_mW, acc_TE, acc_TE_std, '--^b', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'CapSize', 6, 'DisplayName', 'TE Accidentals');
errorbar(powers_TM_mW, rawCC_TM, rawCC_TM_std, '-or', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'CapSize', 6, 'DisplayName', 'TM Raw Coincidences');
errorbar(powers_TM_mW, acc_TM, acc_TM_std, '--^r', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'CapSize', 6, 'DisplayName', 'TM Accidentals');

xlabel('Laser Power (mW)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Counts (cps)', 'FontSize', 12, 'FontWeight', 'bold');
title('Raw Coincidences \& Accidentals vs Power', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 11);

%% --- PLOT 2: True Coincidences vs Power (mW) with Quadratic Fit ---
figure('Name', 'True Coincidences vs Power (mW)', 'Color', 'w');
hold on; box on; grid on;

errorbar(powers_TE_mW, trueCC_TE, trueCC_TE_std, 'ob', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'CapSize', 6, 'LineStyle', 'none', ...
    'DisplayName', 'TE True CC (Raw $-$ Acc)');
errorbar(powers_TM_mW, trueCC_TM, trueCC_TM_std, 'or', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'CapSize', 6, 'LineStyle', 'none', ...
    'DisplayName', 'TM True CC (Raw $-$ Acc)');

fit_eqn     = fittype('a*x^2 + b*x + c');
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, 1, 1]);

if length(powers_TE_mW) >= 3
    [fit_TE, ~] = fit(powers_TE_mW', trueCC_TE', fit_eqn, fit_options);
    x_TE_fit = linspace(min(powers_TE_mW), max(powers_TE_mW), 100);
    plot(x_TE_fit, fit_TE(x_TE_fit), 'b-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('TE Fit: $a=%.3f,\\, b=%.3f,\\, c=%.3f$', ...
        fit_TE.a, fit_TE.b, fit_TE.c));
end

if length(powers_TM_mW) >= 3
    [fit_TM, ~] = fit(powers_TM_mW', trueCC_TM', fit_eqn, fit_options);
    x_TM_fit = linspace(min(powers_TM_mW), max(powers_TM_mW), 100);
    plot(x_TM_fit, fit_TM(x_TM_fit), 'r-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('TM Fit: $a=%.3f,\\, b=%.3f,\\, c=%.3f$', ...
        fit_TM.a, fit_TM.b, fit_TM.c));
end

xlabel('Laser Power (mW)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('True Coincidences (cps)', 'FontSize', 12, 'FontWeight', 'bold');
title('SFWM True Coincidences vs Laser Power', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 11);

%% --- PLOT 3: CAR vs Power (mW) ---
figure('Name', 'CAR vs Power (mW)', 'Color', 'w');
hold on; box on; grid on;

errorbar(powers_TE_mW, car_TE, car_TE_std, '-ob', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'CapSize', 6, 'DisplayName', 'TE Polarization');
errorbar(powers_TM_mW, car_TM, car_TM_std, '-or', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r', 'CapSize', 6, 'DisplayName', 'TM Polarization');

xlabel('Laser Power (mW)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CAR', 'FontSize', 12, 'FontWeight', 'bold');
title('Coincidence-to-Accidental Ratio (CAR) vs Power', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11);

% =========================================================================
%% SECTION 6: PLOTS — dBm POWER AXIS
% =========================================================================

%% --- PLOT 4: Raw Coincidences & Accidentals vs Power (dBm) ---
figure('Name', 'Raw CC and Accidentals vs Power (dBm)', 'Color', 'w');
hold on; box on; grid on;

errorbar(powers_TE_dBm, rawCC_TE, rawCC_TE_std, '-ob', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'CapSize', 6, 'DisplayName', 'TE Raw Coincidences');
errorbar(powers_TE_dBm, acc_TE, acc_TE_std, '--^b', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'CapSize', 6, 'DisplayName', 'TE Accidentals');
errorbar(powers_TM_dBm, rawCC_TM, rawCC_TM_std, '-or', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'CapSize', 6, 'DisplayName', 'TM Raw Coincidences');
errorbar(powers_TM_dBm, acc_TM, acc_TM_std, '--^r', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'CapSize', 6, 'DisplayName', 'TM Accidentals');

xlabel('Laser Power (dBm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Counts (cps)', 'FontSize', 12, 'FontWeight', 'bold');
title('Raw Coincidences \& Accidentals vs Power', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 11);

%% --- PLOT 5: True Coincidences vs Power (dBm) ---
figure('Name', 'True Coincidences vs Power (dBm)', 'Color', 'w');
hold on; box on; grid on;

errorbar(powers_TE_dBm, trueCC_TE, trueCC_TE_std, 'ob', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'CapSize', 6, 'LineStyle', 'none', ...
    'DisplayName', 'TE True CC (Raw $-$ Acc)');
errorbar(powers_TM_dBm, trueCC_TM, trueCC_TM_std, 'or', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'CapSize', 6, 'LineStyle', 'none', ...
    'DisplayName', 'TM True CC (Raw $-$ Acc)');

xlabel('Laser Power (dBm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('True Coincidences (cps)', 'FontSize', 12, 'FontWeight', 'bold');
title('SFWM True Coincidences vs Laser Power', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 11);

%% --- PLOT 6: CAR vs Power (dBm) ---
figure('Name', 'CAR vs Power (dBm)', 'Color', 'w');
hold on; box on; grid on;

errorbar(powers_TE_dBm, car_TE, car_TE_std, '-ob', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'CapSize', 6, 'DisplayName', 'TE Polarization');
errorbar(powers_TM_dBm, car_TM, car_TM_std, '-or', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r', 'CapSize', 6, 'DisplayName', 'TM Polarization');

xlabel('Laser Power (dBm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CAR', 'FontSize', 12, 'FontWeight', 'bold');
title('Coincidence-to-Accidental Ratio (CAR) vs Power', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11);