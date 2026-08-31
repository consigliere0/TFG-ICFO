% =========================================================================
% TFG - CAR / COINCIDENCE RATES: RESULTS PIPELINE WITH ERROR BARS
% =========================================================================
% ÍNDEX:
%   SECTION 1: Settings & color palette
%   SECTION 2: SPAD data loading (full power range, TE/TM)
%   SECTION 3: Error-bar computation (Poisson, per supervisor's formula)
%   SECTION 4: SPAD Figure 1 - Coincidence rates vs power (accid. + true)
%   SECTION 5: SPAD Figure 2 - CAR vs power
%   SECTION 6: SNSPD data loading (decay-tail range, TE/TM; bell excluded)
%   SECTION 7: SNSPD Figure - HH vs VV counts per pumping (axis NOT
%              calibrated - index units only, see caveat below)
%   SECTION 8: SNSPD CAR (optional, gated by ENABLE_SNSPD_CAR flag)
% =========================================================================

%% =========================================================================
% SECTION 1: Settings & color palette
% =========================================================================
clear; close all; format long g
set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

% --- Acquisition parameters ---
T_acq_SPAD  = 60;            % [s] acquisition time per point, SPAD run
tau_SPAD    = 1000 * 1e-12;  % [s] coincidence window, SPAD run
                              % [PLACEHOLDER: confirm 1 ns, not 400 ps -
                              %  differs from your stated Methods value]

T_acq_SNSPD = 60;            % [PLACEHOLDER: confirm acquisition time for
                              %  the entr5 (SNSPD) run - assumed same as SPAD]
tau_SNSPD   = 400 * 1e-12;   % [s] matches Methods (400 ps)

% --- Color palette ---
c_acc_TE  = [0.40, 0.78, 0.50];   % green, TE accidentals
c_acc_TM  = [0.20, 0.55, 0.35];   % dark green, TM accidentals
c_true_TE = [0.35, 0.60, 0.90];   % blue, TE true coincidences
c_true_TM = [0.15, 0.35, 0.65];   % dark blue, TM true coincidences
c_CAR_TE  = [0.25, 0.75, 0.72];   % teal
c_CAR_TM  = [0.97, 0.55, 0.20];   % orange

c_HH = [0.55, 0.30, 0.75];  % purple
c_VV = [0.85, 0.35, 0.35];  % red
c_HV = [0.25, 0.55, 0.75];  % blue
c_VH = [0.80, 0.55, 0.20];  % ochre


%outFolder = '../Tests/timetagger/Results_Figures';
outFolder = 'results_figures';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

%% =========================================================================
% SECTION 2: SPAD data loading (full power range, TE/TM)
% =========================================================================
data_path = '../../Tests/timetagger/wvg2/data5/';
files     = dir(fullfile(data_path, '*2026-*'));

manual_powers_TE_dBm = [-21.84, -18.74, -16.74, -14.84, -12.94, -11.04, -8.94, -6.74, -4.94, -2.84, -0.44, 1.66, 3.66, 5.86, 7.86, 10.16, 11.66];
manual_powers_TM_dBm = [-22.64, -19.84, -17.14, -15.24, -13.04, -11.14, -9.14, -7.14, -5.04, -3.14, -1.04, 0.76, 2.76, 4.56, 6.56, 8.46, 10.46, 11.36];

powers_TE_dBm = []; powers_TM_dBm = [];
rawCC_TE  = []; acc_TE  = []; trueCC_TE = []; car_TE = [];
rawCC_TM  = []; acc_TM  = []; trueCC_TM = []; car_TM = [];

fprintf('Loading SPAD CAR data...\n');
for i = 1:length(files)
    fname          = files(i).name;
    full_file_path = fullfile(files(i).folder, fname);

    token_TE = regexp(fname, 'TE(\d+)', 'tokens');
    token_TM = regexp(fname, 'TM(\d+)', 'tokens');
    if isempty(token_TE) && isempty(token_TM), continue; end

    try
        data = readtable(full_file_path, 'FileType', 'text', 'Delimiter', ',', ...
            'NumHeaderLines', 9, 'VariableNamingRule', 'preserve');
    catch
        warning('Could not read: %s', full_file_path); continue;
    end
    if size(data, 2) < 8
        warning('Too few columns: %s', fname); continue;
    end

    CH1_val   = data{1, 3};
    CH2_val   = data{1, 4};
    Coinc_val = data{1, 8};

    accidentals       = CH1_val * CH2_val * tau_SPAD;
    true_coincidences = Coinc_val - accidentals;
    CAR               = true_coincidences / accidentals;

    if ~isempty(token_TE)
        file_idx  = str2double(token_TE{1}{1});
        powers_TE_dBm(end+1) = manual_powers_TE_dBm(file_idx + 1); %#ok<*SAGROW>
        rawCC_TE(end+1)      = Coinc_val;
        acc_TE(end+1)        = accidentals;
        trueCC_TE(end+1)     = true_coincidences;
        car_TE(end+1)        = CAR;
    elseif ~isempty(token_TM)
        file_idx  = str2double(token_TM{1}{1});
        powers_TM_dBm(end+1) = manual_powers_TM_dBm(file_idx + 1);
        rawCC_TM(end+1)      = Coinc_val;
        acc_TM(end+1)        = accidentals;
        trueCC_TM(end+1)     = true_coincidences;
        car_TM(end+1)        = CAR;
    end
end
fprintf('Loaded: %d TE points, %d TM points\n', length(powers_TE_dBm), length(powers_TM_dBm));

[powers_TE_dBm, sIdx_TE] = sort(powers_TE_dBm);
rawCC_TE = rawCC_TE(sIdx_TE); acc_TE = acc_TE(sIdx_TE);
trueCC_TE = trueCC_TE(sIdx_TE); car_TE = car_TE(sIdx_TE);

[powers_TM_dBm, sIdx_TM] = sort(powers_TM_dBm);
rawCC_TM = rawCC_TM(sIdx_TM); acc_TM = acc_TM(sIdx_TM);
trueCC_TM = trueCC_TM(sIdx_TM); car_TM = car_TM(sIdx_TM);

powers_TE_mW = 10.^(powers_TE_dBm / 10);
powers_TM_mW = 10.^(powers_TM_dBm / 10);

%% =========================================================================
% SECTION 3: Error-bar computation (Poisson, per supervisor's formula)
% =========================================================================
% sigma_R = sqrt(R / T), for a rate R [cps] measured over acquisition
% time T [s]. Applied to the TRUE-coincidence rate only (per your
% supervisor's simplification); accidentals are treated as effectively
% error-free since they are a deterministic product of two well-sampled
% singles rates.

err_trueCC_TE = poissonRateError(trueCC_TE, T_acq_SPAD);
err_trueCC_TM = poissonRateError(trueCC_TM, T_acq_SPAD);

% CAR = trueCC / acc  =>  sigma_CAR = sigma_trueCC / acc
err_car_TE = err_trueCC_TE ./ acc_TE;
err_car_TM = err_trueCC_TM ./ acc_TM;

%% =========================================================================
% SECTION 4: SPAD Figure 1 - Coincidence rates vs power (accid. + true)
% =========================================================================
fig1 = figure('Name', 'SPAD Coincidence Rates vs Power', 'Color', 'w', ...
              'Position', [80 120 800 520]);
hold on; box on; grid on; grid minor;

errorbar(powers_TE_mW, acc_TE, zeros(size(acc_TE)), 'o', ...
    'Color', c_acc_TE, 'MarkerFaceColor', c_acc_TE, 'MarkerSize', 6, ...
    'LineWidth', 1.2, 'DisplayName', 'TE Accidentals');
errorbar(powers_TM_mW, acc_TM, zeros(size(acc_TM)), 's', ...
    'Color', c_acc_TM, 'MarkerFaceColor', c_acc_TM, 'MarkerSize', 6, ...
    'LineWidth', 1.2, 'DisplayName', 'TM Accidentals');

errorbar(powers_TE_mW, trueCC_TE, err_trueCC_TE, 'o', ...
    'Color', c_true_TE, 'MarkerFaceColor', c_true_TE, 'MarkerSize', 6, ...
    'LineWidth', 1.2, 'DisplayName', 'TE True coincidences');
errorbar(powers_TM_mW, trueCC_TM, err_trueCC_TM, 's', ...
    'Color', c_true_TM, 'MarkerFaceColor', c_true_TM, 'MarkerSize', 6, ...
    'LineWidth', 1.2, 'DisplayName', 'TM True coincidences');

set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickDir', 'out');
xlabel('Pump power (mW)', 'FontSize', 13);
ylabel('Count rate (cps)', 'FontSize', 13);
title('SPAD coincidence rates vs. pump power', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2);
hold off;
saveas(fig1, fullfile(outFolder, 'SPAD_Coincidence_Rates.png'));

%% =========================================================================
% SECTION 5: SPAD Figure 2 - CAR vs power
% =========================================================================
fig2 = figure('Name', 'SPAD CAR vs Power', 'Color', 'w', 'Position', [80 100 800 480]);
hold on; box on; grid on; grid minor;

errorbar(powers_TE_mW, car_TE, err_car_TE, 'o-', ...
    'Color', c_CAR_TE, 'MarkerFaceColor', c_CAR_TE, 'MarkerSize', 7, ...
    'LineWidth', 1.5, 'DisplayName', 'TE');
errorbar(powers_TM_mW, car_TM, err_car_TM, 's-', ...
    'Color', c_CAR_TM, 'MarkerFaceColor', c_CAR_TM, 'MarkerSize', 7, ...
    'LineWidth', 1.5, 'DisplayName', 'TM');

set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickDir', 'out');
xlabel('Pump power (mW)', 'FontSize', 13);
ylabel('CAR', 'FontSize', 13);
title('Coincidence-to-accidental ratio vs. pump power (SPAD)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 11);
hold off;
saveas(fig2, fullfile(outFolder, 'SPAD_CAR_vs_Power.png'));

fprintf('\n--- SPAD CAR summary ---\n');
fprintf('TE: max CAR = %.2f +/- %.2f at %.2f dBm\n', max(car_TE), err_car_TE(car_TE==max(car_TE)), powers_TE_dBm(car_TE==max(car_TE)));
fprintf('TM: max CAR = %.2f +/- %.2f at %.2f dBm\n', max(car_TM), err_car_TM(car_TM==max(car_TM)), powers_TM_dBm(car_TM==max(car_TM)));

%% =========================================================================
% SECTION 6: SNSPD data loading (decay-tail range, TE/TM; bell excluded)
% =========================================================================
% NOTE: the 'bell' dataset from your original script is deliberately
% excluded per your instruction (it was not a true Bell-state measurement).

data_path_snspd = '../../Tests/timetagger/wvg2/entr5';

[raw_HH_te, raw_VV_te, raw_HV_te, raw_VH_te, true_HH_te, true_VV_te, true_HV_te, true_VH_te] = ...
    loadSNSPDPolData(data_path_snspd, '*TE*.txt', tau_SNSPD);

[raw_HH_tm, raw_VV_tm, raw_HV_tm, raw_VH_tm, true_HH_tm, true_VV_tm, true_HV_tm, true_VH_tm] = ...
    loadSNSPDPolData(data_path_snspd, '*TM*.txt', tau_SNSPD);

% [PLACEHOLDER: replace with real coupled power once VOA-to-power mapping
% is calibrated, exactly as done for the SPAD manual_powers_* arrays above]
idx_TE = 0:(length(true_HH_te)-1);
idx_TM = 0:(length(true_HH_tm)-1);

err_HH_te = poissonRateError(true_HH_te, T_acq_SNSPD);
err_VV_te = poissonRateError(true_VV_te, T_acq_SNSPD);
err_HV_te = poissonRateError(true_HV_te, T_acq_SNSPD);
err_VH_te = poissonRateError(true_VH_te, T_acq_SNSPD);

err_HH_tm = poissonRateError(true_HH_tm, T_acq_SNSPD);
err_VV_tm = poissonRateError(true_VV_tm, T_acq_SNSPD);
err_HV_tm = poissonRateError(true_HV_tm, T_acq_SNSPD);
err_VH_tm = poissonRateError(true_VH_tm, T_acq_SNSPD);

%% =========================================================================
% SECTION 7: SNSPD Figure - all polarization coincidence rates
% =========================================================================

fig3 = figure('Name', 'SNSPD Coincidence Rates', 'Color', 'w', ...
              'Position', [80 100 1100 450]);

% -------------------------------------------------------------------------
% TE pumping
% -------------------------------------------------------------------------
subplot(1,2,1);
hold on; box on; grid on;

errorbar(idx_TE, true_HH_te, err_HH_te, 'o-', ...
    'Color', c_HH, ...
    'MarkerFaceColor', c_HH, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'HH');

errorbar(idx_TE, true_VV_te, err_VV_te, 's-', ...
    'Color', c_VV, ...
    'MarkerFaceColor', c_VV, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'VV');

errorbar(idx_TE, true_HV_te, err_HV_te, '^-', ...
    'Color', c_HV, ...
    'MarkerFaceColor', c_HV, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'HV');

errorbar(idx_TE, true_VH_te, err_VH_te, 'd-', ...
    'Color', c_VH, ...
    'MarkerFaceColor', c_VH, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'VH');

xlabel('Measurement index', 'FontSize', 11);
ylabel('True coincidence rate (cps)', 'FontSize', 11);
title('TE pumping', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
set(gca, 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 0.8);

% -------------------------------------------------------------------------
% TM pumping
% -------------------------------------------------------------------------
subplot(1,2,2);
hold on; box on; grid on;

errorbar(idx_TM, true_HH_tm, err_HH_tm, 'o-', ...
    'Color', c_HH, ...
    'MarkerFaceColor', c_HH, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'HH');

errorbar(idx_TM, true_VV_tm, err_VV_tm, 's-', ...
    'Color', c_VV, ...
    'MarkerFaceColor', c_VV, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'VV');

errorbar(idx_TM, true_HV_tm, err_HV_tm, '^-', ...
    'Color', c_HV, ...
    'MarkerFaceColor', c_HV, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'HV');

errorbar(idx_TM, true_VH_tm, err_VH_tm, 'd-', ...
    'Color', c_VH, ...
    'MarkerFaceColor', c_VH, ...
    'LineWidth', 1.3, ...
    'DisplayName', 'VH');

xlabel('Measurement index', 'FontSize', 11);
ylabel('True coincidence rate (cps)', 'FontSize', 11);
title('TM pumping', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
set(gca, 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 0.8);

sgtitle('SNSPD coincidence rates', ...
        'FontSize', 13, 'Interpreter', 'latex');

hold off;

saveas(fig3, fullfile(outFolder, 'SNSPD_Coincidence_Rates.png'));
fprintf('\n--- SNSPD diagnostic: which channel dominates? ---\n');
fprintf('TE pumping: mean HH = %.1f cps, mean VV = %.1f cps  (expect HH >> VV if TE=H convention holds)\n', mean(true_HH_te), mean(true_VV_te));
fprintf('TM pumping: mean HH = %.1f cps, mean VV = %.1f cps  (expect VV >> HH if TM=V convention holds)\n', mean(true_HH_tm), mean(true_VV_tm));

%% =========================================================================
% SECTION 8: SNSPD CAR (optional, gated - do not enable until axis fixed)
% =========================================================================
ENABLE_SNSPD_CAR = false;  % [PLACEHOLDER: set true only once power axis is calibrated]

if ENABLE_SNSPD_CAR
    % Uses the dominant channel per pumping, as determined by Section 7's
    % diagnostic printout above - EDIT these two lines once confirmed:
    dominant_TE = true_HH_te; acc_dominant_TE = raw_HH_te - true_HH_te; %#ok<UNRCH>
    dominant_TM = true_VV_tm; acc_dominant_TM = raw_VV_tm - true_VV_tm;

    car_snspd_TE = dominant_TE ./ acc_dominant_TE;
    car_snspd_TM = dominant_TM ./ acc_dominant_TM;
    err_car_snspd_TE = poissonRateError(dominant_TE, T_acq_SNSPD) ./ acc_dominant_TE;
    err_car_snspd_TM = poissonRateError(dominant_TM, T_acq_SNSPD) ./ acc_dominant_TM;

    fig4 = figure('Name', 'SNSPD CAR (decay tail)', 'Color', 'w', 'Position', [80 100 800 480]);
    hold on; box on; grid on; grid minor;
    errorbar(idx_TE, car_snspd_TE, err_car_snspd_TE, 'o-', 'Color', c_CAR_TE, ...
        'MarkerFaceColor', c_CAR_TE, 'LineWidth', 1.5, 'DisplayName', 'TE');
    errorbar(idx_TM, car_snspd_TM, err_car_snspd_TM, 's-', 'Color', c_CAR_TM, ...
        'MarkerFaceColor', c_CAR_TM, 'LineWidth', 1.5, 'DisplayName', 'TM');
    xlabel('Measurement index (power axis NOT calibrated)', 'FontSize', 12);
    ylabel('CAR', 'FontSize', 13);
    title('SNSPD CAR, decay-tail regime -- UNCALIBRATED AXIS', 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 11);
    hold off;
    saveas(fig4, fullfile(outFolder, 'SNSPD_CAR_UNCALIBRATED.png'));
end

fprintf('\nAll figures saved to: %s\n', outFolder);


%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function sigma_R = poissonRateError(R, T)
% Poisson standard error on a rate R [cps] measured over acquisition
% time T [s]: recover raw count N = R*T, sigma_N = sqrt(N), then
% sigma_R = sigma_N / T = sqrt(R/T). Negative/zero rates guarded to 0.
    R = max(R, 0);
    sigma_R = sqrt(R ./ T);
end


function [raw_HH, raw_VV, raw_HV, raw_VH, true_HH, true_VV, true_HV, true_VH] = ...
    loadSNSPDPolData(data_path, filePattern, tau)
% Loads and processes one polarization run (TE or TM) of four-channel
% SNSPD coincidence data, following the original column mapping:
%   col 3=HA singles, col4=VA singles, col5=HB singles, col6=VB singles
%   col 8=HH raw CC, col9=VV raw CC, col10=HV raw CC, col11=VH raw CC

    files = dir(fullfile(data_path, filePattern));
    raw_HH = []; raw_VV = []; raw_HV = []; raw_VH = [];
    ch1 = []; ch2 = []; ch3 = []; ch4 = [];

    for k = 1:length(files)
        file_path = fullfile(files(k).folder, files(k).name);
        data = readmatrix(file_path, 'NumHeaderLines', 12);

        ch1(end+1) = data(1,3); %#ok<AGROW> % HA
        ch2(end+1) = data(1,4);              % VA
        ch3(end+1) = data(1,5);              % HB
        ch4(end+1) = data(1,6);              % VB

        raw_HH(end+1) = data(1,8);
        raw_VV(end+1) = data(1,9);
        raw_HV(end+1) = data(1,10);
        raw_VH(end+1) = data(1,11);
    end

    acc_HH = ch1 .* ch3 * tau;
    acc_HV = ch1 .* ch4 * tau;
    acc_VV = ch2 .* ch4 * tau;
    acc_VH = ch2 .* ch3 * tau;

    true_HH = raw_HH - acc_HH;
    true_HV = raw_HV - acc_HV;
    true_VV = raw_VV - acc_VV;
    true_VH = raw_VH - acc_VH;
end