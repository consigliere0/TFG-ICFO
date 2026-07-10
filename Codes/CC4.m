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
data_path = '../Tests/timetagger/wvg2/data5/';
files     = dir(fullfile(data_path, '*2026-*'));

manual_powers_TE_dBm = [-21.84, -18.74, -16.74, -14.84, -12.94, -11.04, -8.94, -6.74, -4.94, -2.84, -0.44, 1.66, 3.66, 5.86, 7.86, 10.16, 11.66];
manual_powers_TM_dBm = [-22.64, -19.84, -17.14, -15.24, -13.04, -11.14, -9.14, -7.14, -5.04, -3.14, -1.04, 0.76, 2.76, 4.56, 6.56, 8.46, 10.46, 11.36];

powers_TE_dBm = []; powers_TM_dBm = [];

rawCC_TE  = []; acc_TE  = []; trueCC_TE = []; car_TE = [];
rawCC_TM  = []; acc_TM  = []; trueCC_TM = []; car_TM = [];

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
    % Single datapoint per file — read the one value directly, no mean/std needed
    CH1_val   = data{1, 3};
    CH2_val   = data{1, 4};
    Coinc_val = data{1, 8};

    accidentals       = CH1_val * CH2_val * tau;
    true_coincidences = Coinc_val - accidentals;
    CAR               = true_coincidences / accidentals;

    if ~isempty(token_TE)
        file_idx  = str2double(token_TE{1}{1});
        power_dBm = manual_powers_TE_dBm(file_idx + 1);
        powers_TE_dBm(end+1) = power_dBm;
        rawCC_TE(end+1)      = Coinc_val;
        acc_TE(end+1)        = accidentals;
        trueCC_TE(end+1)     = true_coincidences;
        car_TE(end+1)        = CAR;
        fprintf('TE: %s | Idx %d | %.2f dBm | Raw=%.2f | Acc=%.4f | True=%.2f | CAR=%.2f\n', ...
            fname, file_idx, power_dBm, Coinc_val, accidentals, true_coincidences, CAR);

    elseif ~isempty(token_TM)
        file_idx  = str2double(token_TM{1}{1});
        power_dBm = manual_powers_TM_dBm(file_idx + 1);
        powers_TM_dBm(end+1) = power_dBm;
        rawCC_TM(end+1)      = Coinc_val;
        acc_TM(end+1)        = accidentals;
        trueCC_TM(end+1)     = true_coincidences;
        car_TM(end+1)        = CAR;
        fprintf('TM: %s | Idx %d | %.2f dBm | Raw=%.2f | Acc=%.4f | True=%.2f | CAR=%.2f\n', ...
            fname, file_idx, power_dBm, Coinc_val, accidentals, true_coincidences, CAR);
    end
end

fprintf('\nFiles parsed: %d TE, %d TM\n', length(powers_TE_dBm), length(powers_TM_dBm));

% =========================================================================
%% SECTION 4: SORT
% =========================================================================
[powers_TE_dBm, sIdx_TE] = sort(powers_TE_dBm);
rawCC_TE  = rawCC_TE(sIdx_TE);
acc_TE    = acc_TE(sIdx_TE);
trueCC_TE = trueCC_TE(sIdx_TE);
car_TE    = car_TE(sIdx_TE);

[powers_TM_dBm, sIdx_TM] = sort(powers_TM_dBm);
rawCC_TM  = rawCC_TM(sIdx_TM);
acc_TM    = acc_TM(sIdx_TM);
trueCC_TM = trueCC_TM(sIdx_TM);
car_TM    = car_TM(sIdx_TM);

powers_TE_mW = 10.^(powers_TE_dBm / 10);
powers_TM_mW = 10.^(powers_TM_dBm / 10);

disp('Sorted. Generating plots...');

% =========================================================================
%% SECTION 5: COLOUR PALETTE
% =========================================================================
c_raw_data  = [0.92, 0.45, 0.45];
c_raw_fit   = [0.97, 0.75, 0.75];

c_acc_data  = [0.40, 0.78, 0.50];
c_acc_fit   = [0.75, 0.93, 0.80];

c_true_data = [0.35, 0.60, 0.90];
c_true_fit  = [0.72, 0.85, 0.97];

c_TE_CAR    = [0.25, 0.75, 0.72];
c_TM_CAR    = [0.97, 0.72, 0.25];

% =========================================================================
%% SECTION 6: PLOT 1 — Coincidences: two subplots side by side (dBm | mW)
% =========================================================================
fit_eqn     = fittype('a*x^2 + b*x + c');
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, 1, 1]);

figure('Name', 'All Coincidences vs Power', 'Color', 'w', 'Position', [80 120 1400 520]);

for subplot_idx = 1:2
    if subplot_idx == 1
        ax    = subplot(1, 2, 1);
        px_TE = powers_TE_dBm;
        px_TM = powers_TM_dBm;
        xlab  = 'Pump Power (dBm)';
        tstr  = 'Coincidence Rates vs Pump Power (dBm)';
    else
        ax    = subplot(1, 2, 2);
        px_TE = powers_TE_mW;
        px_TM = powers_TM_mW;
        xlab  = 'Pump Power (mW)';
        tstr  = 'Coincidence Rates vs Pump Power (mW)';
    end

    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); grid(ax, 'minor');

    % % ---- RAW CC ----
    % plot(ax, px_TE, rawCC_TE, 'o', ...
    %     'Color', c_raw_data, 'MarkerFaceColor', c_raw_data, 'MarkerSize', 7, ...
    %     'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', 'TE Raw CC');
    % plot(ax, px_TM, rawCC_TM, 's', ...
    %     'Color', c_raw_data, 'MarkerFaceColor', c_raw_data, 'MarkerSize', 7, ...
    %     'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', 'TM Raw CC');

    % ---- ACCIDENTALS ----
    plot(ax, px_TE, acc_TE, 'o', ...
        'Color', c_acc_data, 'MarkerFaceColor', c_acc_data, 'MarkerSize', 7, ...
        'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', 'TE Accidentals');
    plot(ax, px_TM, acc_TM, 's', ...
        'Color', c_acc_data, 'MarkerFaceColor', c_acc_data, 'MarkerSize', 7, ...
        'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', 'TM Accidentals');

    % if length(px_TE) >= 3
    %     [f, ~] = fit(px_TE', acc_TE', fit_eqn, fit_options);
    %     xf = linspace(min(px_TE), max(px_TE), 200);
    %     plot(ax, xf, f(xf), '-', 'Color', c_acc_fit, 'LineWidth', 2.0, ...
    %         'DisplayName', 'TE Acc fit');
    % end
    % if length(px_TM) >= 3
    %     [f, ~] = fit(px_TM', acc_TM', fit_eqn, fit_options);
    %     xf = linspace(min(px_TM), max(px_TM), 200);
    %     plot(ax, xf, f(xf), '--', 'Color', c_acc_fit, 'LineWidth', 2.0, ...
    %         'DisplayName', 'TM Acc fit');
    % end

    % ---- TRUE CC ----
    plot(ax, px_TE, trueCC_TE, 'o', ...
        'Color', c_true_data, 'MarkerFaceColor', c_true_data, 'MarkerSize', 7, ...
        'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', 'TE True CC');
    plot(ax, px_TM, trueCC_TM, 's', ...
        'Color', [0.92, 0.45, 0.45], 'MarkerFaceColor', [0.92, 0.45, 0.45], 'MarkerSize', 7, ...
        'LineWidth', 1.2, 'LineStyle', 'none', 'DisplayName', 'TM True CC');

    xlabel(ax, xlab,               'FontSize', 13);
    ylabel(ax, 'Count Rate (cps)', 'FontSize', 13);
    title(ax,  tstr,               'FontSize', 13);
    legend(ax, 'Location', 'northwest', 'FontSize', 9, 'NumColumns', 2);
    set(ax, 'FontSize', 11, 'LineWidth', 0.8);
end

sgtitle('Coincidence Rates vs Pump Power', 'FontSize', 14, 'Interpreter', 'latex');

% =========================================================================
%% SECTION 7: PLOT 2 — CAR: two subplots side by side (dBm | mW)
% =========================================================================
figure('Name', 'CAR vs Power', 'Color', 'w', 'Position', [80 100 1000 430]);

% ---- Left subplot: dBm ----
% ax1 = subplot(1, 2, 1);
% hold on; box on; grid on; grid minor;
% 
% plot(powers_TE_dBm, car_TE, 'o-', ...
%     'Color', c_TE_CAR, 'MarkerFaceColor', c_TE_CAR, 'MarkerSize', 8, ...
%     'LineWidth', 1.5, 'DisplayName', 'TE Polarization');
% plot(powers_TM_dBm, car_TM, 's-', ...
%     'Color', c_TM_CAR, 'MarkerFaceColor', c_TM_CAR, 'MarkerSize', 8, ...
%     'LineWidth', 1.5, 'DisplayName', 'TM Polarization');
% 
% xlabel('Pump Power (dBm)', 'FontSize', 13);
% ylabel('CAR',              'FontSize', 13);
% title('CAR vs Power (dBm)', 'FontSize', 13);
% legend('Location', 'northeast', 'FontSize', 11);
% set(ax1, 'FontSize', 11, 'LineWidth', 0.8);

% ---- Right subplot: mW ----
%ax2 = subplot(1, 2, 2);
hold on; box on; grid on; grid minor;

plot(powers_TE_mW, car_TE, 'o-', ...
    'Color', c_TE_CAR, 'MarkerFaceColor', c_TE_CAR, 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'DisplayName', 'TE Polarization');
plot(powers_TM_mW, car_TM, 's-', ...
    'Color', c_TM_CAR, 'MarkerFaceColor', c_TM_CAR, 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'DisplayName', 'TM Polarization');

xlabel('Pump Power (mW)', 'FontSize', 13);
ylabel('CAR',             'FontSize', 13);
title('CAR vs Power (mW)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 11);
%set(ax2, 'FontSize', 11, 'LineWidth', 0.8);

%sgtitle('Coincidence-to-Accidental Ratio vs Pump Power', 'FontSize', 14, 'Interpreter', 'latex');