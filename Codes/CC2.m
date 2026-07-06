% =========================================================================
%% TFG - COINCIDENCE COUNTS & CAR DATA PROCESSING
% =========================================================================
%%                     SECTION 1: INITIAL SETTINGS
% =========================================================================
clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

tau = 1000 * 1e-12;

% =========================================================================
%% SECTION 2: LOAD DATA & MANUAL POWER MAPPING
% =========================================================================
data_path = '../Tests/timetagger/wvg2/data3/';
% 5 repeats
files     = dir(fullfile(data_path, '*2026-*'));

manual_powers_TE_dBm = [ 13.38, 9.58, -0.42, -10.42, -13.72, 4.58 ];

manual_powers_TM_dBm = [ -5, -10, -15, -20, -25, -30 ];

powers_TE_dBm = []; powers_TM_dBm = [];

rawCC_TE = []; rawCC_TE_std = []; acc_TE = []; acc_TE_std = [];
trueCC_TE = []; trueCC_TE_std = []; car_TE = []; car_TE_std = [];

rawCC_TM = []; rawCC_TM_std = []; acc_TM = []; acc_TM_std = [];
trueCC_TM = []; trueCC_TM_std = []; car_TM = []; car_TM_std = [];

disp('Processing data files...');

for i = 1:length(files)
    fname          = files(i).name;
    full_file_path = fullfile(files(i).folder, fname);

    token_TE = regexp(fname, 'TE(\d+)', 'tokens');
    token_TM = regexp(fname, 'TM(\d+)', 'tokens');

    if isempty(token_TE) && isempty(token_TM)
        fprintf('Skipping: %s\n', fname); continue;
    end

    try
        data = readtable(full_file_path, 'FileType', 'text', 'Delimiter', ',', ...
            'NumHeaderLines', 9, 'VariableNamingRule', 'preserve');
    catch
        warning('Could not read: %s', full_file_path); continue;
    end

    if size(data, 2) < 8
        warning('Too few columns: %s', fname); continue;
    end

    % =====================================================================
    %% SECTION 3: MAGNITUDE CALCULATION
    % =====================================================================
    CH1_vec   = data{:, 3};
    CH2_vec   = data{:, 4};
    Coinc_vec = data{:, 8};

    CH1_mean       = mean(CH1_vec,   'omitnan');
    CH2_mean       = mean(CH2_vec,   'omitnan');
    Coinc_raw_mean = mean(Coinc_vec, 'omitnan');

    CH1_std       = std(CH1_vec,   0, 'omitnan');
    CH2_std       = std(CH2_vec,   0, 'omitnan');
    Coinc_raw_std = std(Coinc_vec, 0, 'omitnan');

    accidentals     = CH1_mean * CH2_mean * tau;
    accidentals_std = sqrt((CH2_mean*tau*CH1_std)^2 + (CH1_mean*tau*CH2_std)^2);

    true_coincidences     = Coinc_raw_mean - accidentals;
    true_coincidences_std = sqrt(Coinc_raw_std^2 + accidentals_std^2);

    CAR     = true_coincidences / accidentals;
    CAR_std = abs(CAR) * sqrt((true_coincidences_std/true_coincidences)^2 + ...
                               (accidentals_std/accidentals)^2);

    if ~isempty(token_TE)
        file_idx  = str2double(token_TE{1}{1});
        power_dBm = manual_powers_TE_dBm(file_idx + 1);
        powers_TE_dBm(end+1)    = power_dBm;
        rawCC_TE(end+1)         = Coinc_raw_mean; rawCC_TE_std(end+1)   = Coinc_raw_std;
        acc_TE(end+1)           = accidentals;    acc_TE_std(end+1)     = accidentals_std;
        trueCC_TE(end+1)        = true_coincidences; trueCC_TE_std(end+1) = true_coincidences_std;
        car_TE(end+1)           = CAR;            car_TE_std(end+1)     = CAR_std;
        fprintf('TE: %s | Idx %d | %.2f dBm | Raw=%.2f | Acc=%.4f | True=%.2f | CAR=%.2f\n', ...
            fname, file_idx, power_dBm, Coinc_raw_mean, accidentals, true_coincidences, CAR);

    elseif ~isempty(token_TM)
        file_idx  = str2double(token_TM{1}{1});
        power_dBm = manual_powers_TM_dBm(file_idx + 1);
        powers_TM_dBm(end+1)    = power_dBm;
        rawCC_TM(end+1)         = Coinc_raw_mean; rawCC_TM_std(end+1)   = Coinc_raw_std;
        acc_TM(end+1)           = accidentals;    acc_TM_std(end+1)     = accidentals_std;
        trueCC_TM(end+1)        = true_coincidences; trueCC_TM_std(end+1) = true_coincidences_std;
        car_TM(end+1)           = CAR;            car_TM_std(end+1)     = CAR_std;
        fprintf('TM: %s | Idx %d | %.2f dBm | Raw=%.2f | Acc=%.4f | True=%.2f | CAR=%.2f\n', ...
            fname, file_idx, power_dBm, Coinc_raw_mean, accidentals, true_coincidences, CAR);
    end
end

fprintf('\nFiles parsed: %d TE, %d TM\n', length(powers_TE_dBm), length(powers_TM_dBm));

% =========================================================================
%% SECTION 4: SORT
% =========================================================================
[powers_TE_dBm, sIdx_TE] = sort(powers_TE_dBm);
rawCC_TE  = rawCC_TE(sIdx_TE);  rawCC_TE_std  = rawCC_TE_std(sIdx_TE);
acc_TE    = acc_TE(sIdx_TE);    acc_TE_std    = acc_TE_std(sIdx_TE);
trueCC_TE = trueCC_TE(sIdx_TE); trueCC_TE_std = trueCC_TE_std(sIdx_TE);
car_TE    = car_TE(sIdx_TE);    car_TE_std    = car_TE_std(sIdx_TE);

[powers_TM_dBm, sIdx_TM] = sort(powers_TM_dBm);
rawCC_TM  = rawCC_TM(sIdx_TM);  rawCC_TM_std  = rawCC_TM_std(sIdx_TM);
acc_TM    = acc_TM(sIdx_TM);    acc_TM_std    = acc_TM_std(sIdx_TM);
trueCC_TM = trueCC_TM(sIdx_TM); trueCC_TM_std = trueCC_TM_std(sIdx_TM);
car_TM    = car_TM(sIdx_TM);    car_TM_std    = car_TM_std(sIdx_TM);

% dBm -> mW for the CAR subplot
powers_TE_mW = 10.^(powers_TE_dBm / 10);
powers_TM_mW = 10.^(powers_TM_dBm / 10);

disp('Sorted. Generating plots...');

% =========================================================================
%% SECTION 5: COLOUR PALETTE
% Each physical quantity gets its own distinct pastel hue.
% Data markers are saturated; fit lines are a lighter tint of the same hue.
% TE and TM share the same hue per quantity, distinguished by marker shape.
% =========================================================================
% Raw CC     — pastel red
c_raw_data  = [0.92, 0.45, 0.45];
c_raw_fit   = [0.97, 0.75, 0.75];

% Accidentals — pastel green
c_acc_data  = [0.40, 0.78, 0.50];
c_acc_fit   = [0.75, 0.93, 0.80];

% True CC    — pastel blue
c_true_data = [0.35, 0.60, 0.90];
c_true_fit  = [0.72, 0.85, 0.97];

% CAR        — TE: pastel teal | TM: pastel amber
c_TE_CAR    = [0.25, 0.75, 0.72];
c_TM_CAR    = [0.97, 0.72, 0.25];

% =========================================================================
%% SECTION 6: PLOT 1 — Superposed Coincidences vs Power (dBm)
%  All six series (TE+TM) × (raw, acc, true) with individual fits
% =========================================================================
fit_eqn     = fittype('a*x^2 + b*x + c');
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, 1, 1]);

figure('Name', 'All Coincidences vs Power', 'Color', 'w', 'Position', [80 120 820 520]);
hold on; box on; grid on; grid minor;

% ---- RAW CC ----
errorbar(powers_TE_dBm, rawCC_TE, rawCC_TE_std, 'o', ...
    'Color', c_raw_data, 'MarkerFaceColor', c_raw_data, 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none', ...
    'DisplayName', 'TE Raw CC');
errorbar(powers_TM_dBm, rawCC_TM, rawCC_TM_std, 's', ...
    'Color', c_raw_data, 'MarkerFaceColor', c_raw_data, 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none', ...
    'DisplayName', 'TM Raw CC');

% if length(powers_TE_dBm) >= 3
%     [f, ~] = fit(powers_TE_dBm', rawCC_TE', fit_eqn, fit_options);
%     xf = linspace(min(powers_TE_dBm), max(powers_TE_dBm), 200);
%     plot(xf, f(xf), '-', 'Color', c_raw_fit, 'LineWidth', 2.0, ...
%         'DisplayName', 'TE Raw fit');
% end
% if length(powers_TM_dBm) >= 3
%     [f, ~] = fit(powers_TM_dBm', rawCC_TM', fit_eqn, fit_options);
%     xf = linspace(min(powers_TM_dBm), max(powers_TM_dBm), 200);
%     plot(xf, f(xf), '--', 'Color', c_raw_fit, 'LineWidth', 2.0, ...
%         'DisplayName', 'TM Raw fit');
% end

% ---- ACCIDENTALS ----
errorbar(powers_TE_dBm, acc_TE, acc_TE_std, 'o', ...
    'Color', c_acc_data, 'MarkerFaceColor', c_acc_data, 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none', ...
    'DisplayName', 'TE Accidentals');
errorbar(powers_TM_dBm, acc_TM, acc_TM_std, 's', ...
    'Color', c_acc_data, 'MarkerFaceColor', c_acc_data, 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none', ...
    'DisplayName', 'TM Accidentals');

if length(powers_TE_dBm) >= 3
    [f, ~] = fit(powers_TE_dBm', acc_TE', fit_eqn, fit_options);
    xf = linspace(min(powers_TE_dBm), max(powers_TE_dBm), 200);
    plot(xf, f(xf), '-', 'Color', c_acc_fit, 'LineWidth', 2.0, ...
        'DisplayName', 'TE Acc fit');
end
if length(powers_TM_dBm) >= 3
    [f, ~] = fit(powers_TM_dBm', acc_TM', fit_eqn, fit_options);
    xf = linspace(min(powers_TM_dBm), max(powers_TM_dBm), 200);
    plot(xf, f(xf), '--', 'Color', c_acc_fit, 'LineWidth', 2.0, ...
        'DisplayName', 'TM Acc fit');
end

% ---- TRUE CC ----
errorbar(powers_TE_dBm, trueCC_TE, trueCC_TE_std, 'o', ...
    'Color', c_true_data, 'MarkerFaceColor', c_true_data, 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none', ...
    'DisplayName', 'TE True CC');
errorbar(powers_TM_dBm, trueCC_TM, trueCC_TM_std, 's', ...
    'Color', c_true_data, 'MarkerFaceColor', c_true_data, 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'CapSize', 5, 'LineStyle', 'none', ...
    'DisplayName', 'TM True CC');

% if length(powers_TE_dBm) >= 3
%     [f, ~] = fit(powers_TE_dBm', trueCC_TE', fit_eqn, fit_options);
%     xf = linspace(min(powers_TE_dBm), max(powers_TE_dBm), 200);
%     plot(xf, f(xf), '-', 'Color', c_true_fit, 'LineWidth', 2.0, ...
%         'DisplayName', 'TE True fit');
% end
% if length(powers_TM_dBm) >= 3
%     [f, ~] = fit(powers_TM_dBm', trueCC_TM', fit_eqn, fit_options);
%     xf = linspace(min(powers_TM_dBm), max(powers_TM_dBm), 200);
%     plot(xf, f(xf), '--', 'Color', c_true_fit, 'LineWidth', 2.0, ...
%         'DisplayName', 'TM True fit');
% end

xlabel('Pump Power (dBm)', 'FontSize', 13);
ylabel('Count Rate (cps)',  'FontSize', 13);
title('Coincidence Rates vs Pump Power', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 9, 'NumColumns', 2);
set(gca, 'FontSize', 11, 'LineWidth', 0.8);

% =========================================================================
%% SECTION 7: PLOT 2 — CAR: two subplots side by side (dBm | mW)
% =========================================================================
figure('Name', 'CAR vs Power', 'Color', 'w', 'Position', [80 100 1000 430]);

% ---- Left subplot: dBm ----
ax1 = subplot(1, 2, 1);
hold on; box on; grid on; grid minor;

errorbar(powers_TE_dBm, car_TE, car_TE_std, 'o-', ...
    'Color', c_TE_CAR, 'MarkerFaceColor', c_TE_CAR, 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'CapSize', 5, 'DisplayName', 'TE Polarization');
errorbar(powers_TM_dBm, car_TM, car_TM_std, 's-', ...
    'Color', c_TM_CAR, 'MarkerFaceColor', c_TM_CAR, 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'CapSize', 5, 'DisplayName', 'TM Polarization');

xlabel('Pump Power (dBm)', 'FontSize', 13);
ylabel('CAR',              'FontSize', 13);
title('CAR vs Power (dBm)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 11);
set(ax1, 'FontSize', 11, 'LineWidth', 0.8);

% ---- Right subplot: mW ----
ax2 = subplot(1, 2, 2);
hold on; box on; grid on; grid minor;

errorbar(powers_TE_mW, car_TE, car_TE_std, 'o-', ...
    'Color', c_TE_CAR, 'MarkerFaceColor', c_TE_CAR, 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'CapSize', 5, 'DisplayName', 'TE Polarization');
errorbar(powers_TM_mW, car_TM, car_TM_std, 's-', ...
    'Color', c_TM_CAR, 'MarkerFaceColor', c_TM_CAR, 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'CapSize', 5, 'DisplayName', 'TM Polarization');

xlabel('Pump Power (mW)', 'FontSize', 13);
ylabel('CAR',             'FontSize', 13);
title('CAR vs Power (mW)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 11);
set(ax2, 'FontSize', 11, 'LineWidth', 0.8);

% Shared supertitle
sgtitle('Coincidence-to-Accidental Ratio vs Pump Power', ...
    'FontSize', 14, 'Interpreter', 'latex');