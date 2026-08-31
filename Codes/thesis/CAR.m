% =========================================================================
% TFG - SPAD CAR / COINCIDENCE RATES / BRIGHTNESS: RESULTS PIPELINE
% =========================================================================
% ÍNDEX:
%   SECTION 1: Settings & color palette
%   SECTION 2: SPAD data loading (full power range, TE/TM)
%   SECTION 3: Poisson error-bar computation
%   SECTION 4: Figure 1 - Coincidence rates vs power (+ quadratic fits)
%   SECTION 5: Figure 2 - CAR vs power (mW)
%   SECTION 6: Figure 3 - CAR vs coincidence rate (true CC, cps)
%   SECTION 7: Brightness - power-law fit to true coincidence rate vs power
% =========================================================================

%% =========================================================================
% SECTION 1: Settings & color palette
% =========================================================================
clear; close all; format long g
set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

T_acq_SPAD = 60;             % [s] acquisition time per point
tau_SPAD   = 1000 * 1e-12;   % [s] coincidence window
                              % [PLACEHOLDER: confirm 1 ns for this run]

% --- Base colors: TE = deep blue, TM = maroon ---
c_TE = [0.10, 0.20, 0.55];
c_TM = [0.55, 0.10, 0.15];

% Lighter tints of the same hue, used for accidentals (keeps TE/TM
% families visually grouped: light = accidentals, dark = true/CAR)
c_TE_light = lightenColor(c_TE, 0.65);
c_TM_light = lightenColor(c_TM, 0.65);

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
% SECTION 3: Poisson error-bar computation
% =========================================================================
% sigma_R = sqrt(R / T) applied to the true-coincidence rate; propagated
% to CAR as sigma_CAR = sigma_trueCC / accidentals (accidentals treated
% as effectively error-free, per your supervisor's simplification).

err_trueCC_TE = poissonRateError(trueCC_TE, T_acq_SPAD);
err_trueCC_TM = poissonRateError(trueCC_TM, T_acq_SPAD);

err_car_TE = err_trueCC_TE ./ acc_TE;
err_car_TM = err_trueCC_TM ./ acc_TM;
% 
% %% =========================================================================
% % SECTION 4: Figure 1 - Coincidence rates vs power (+ quadratic fits)
% % =========================================================================
% fig1 = figure('Name', 'SPAD Coincidence Rates vs Power', 'Color', 'w', ...
%               'Position', [80 120 800 520]);
% hold on; box on; grid on; grid minor;
% 
% % --- Experimental points (in legend) ---
% errorbar(powers_TE_mW, acc_TE, zeros(size(acc_TE)), 'o', ...
%     'Color', c_TE_light, 'MarkerFaceColor', c_TE_light, 'MarkerSize', 6, ...
%     'LineWidth', 1.2, 'DisplayName', 'TE Accidentals');
% errorbar(powers_TM_mW, acc_TM, zeros(size(acc_TM)), 's', ...
%     'Color', c_TM_light, 'MarkerFaceColor', c_TM_light, 'MarkerSize', 6, ...
%     'LineWidth', 1.2, 'DisplayName', 'TM Accidentals');
% errorbar(powers_TE_mW, trueCC_TE, err_trueCC_TE, 'o', ...
%     'Color', c_TE, 'MarkerFaceColor', c_TE, 'MarkerSize', 6, ...
%     'LineWidth', 1.2, 'DisplayName', 'TE True coincidences');
% errorbar(powers_TM_mW, trueCC_TM, err_trueCC_TM, 's', ...
%     'Color', c_TM, 'MarkerFaceColor', c_TM, 'MarkerSize', 6, ...
%     'LineWidth', 1.2, 'DisplayName', 'TM True coincidences');
% 
% % % --- Quadratic fits, transparent, NOT in legend ---
% % plotTransparentQuadFit(powers_TE_mW, acc_TE,    c_TE_light);
% % plotTransparentQuadFit(powers_TM_mW, acc_TM,    c_TM_light);
% % plotTransparentQuadFit(powers_TE_mW, trueCC_TE, c_TE);
% % plotTransparentQuadFit(powers_TM_mW, trueCC_TM, c_TM);
% 
% set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickDir', 'out');
% xlabel('Pump power (mW)', 'FontSize', 13);
% ylabel('Coincidence rate (cps)', 'FontSize', 13);
% legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2);
% hold off;
% saveas(fig1, fullfile(outFolder, 'SPAD_Coincidence_Rates.png'));

%% =========================================================================
% SECTION 4b: Figure 1b - Coincidence rates vs power in dBm
% =========================================================================
fig1_dBm = figure('Name', 'SPAD Coincidence Rates vs Pump Power in dBm', ...
                  'Color', 'w', ...
                  'Position', [120 100 800 520]);
hold on; box on; grid on; grid minor;

% --- Experimental points ---
errorbar(powers_TE_dBm, acc_TE, zeros(size(acc_TE)), 'o', ...
    'Color', c_TE_light, ...
    'MarkerFaceColor', c_TE_light, ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'TE Accidentals');

errorbar(powers_TM_dBm, acc_TM, zeros(size(acc_TM)), 's', ...
    'Color', c_TM_light, ...
    'MarkerFaceColor', c_TM_light, ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'TM Accidentals');

errorbar(powers_TE_dBm, trueCC_TE, err_trueCC_TE, 'o', ...
    'Color', c_TE, ...
    'MarkerFaceColor', c_TE, ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'TE True coincidences');

errorbar(powers_TM_dBm, trueCC_TM, err_trueCC_TM, 's', ...
    'Color', c_TM, ...
    'MarkerFaceColor', c_TM, ...
    'MarkerSize', 6, ...
    'LineWidth', 1.2, ...
    'DisplayName', 'TM True coincidences');

% Optional quadratic fits:
% plotTransparentQuadFit(powers_TE_dBm, acc_TE,    c_TE_light);
% plotTransparentQuadFit(powers_TM_dBm, acc_TM,    c_TM_light);
% plotTransparentQuadFit(powers_TE_dBm, trueCC_TE, c_TE);
% plotTransparentQuadFit(powers_TM_dBm, trueCC_TM, c_TM);

set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickDir', 'out');
xlabel('Pump power (dBm)', 'FontSize', 13);
ylabel('Count rate (cps)', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2);

hold off;
saveas(fig1_dBm, fullfile(outFolder, ...
    'SPAD_Coincidence_Rates_vs_Power_dBm.png'));


%% =========================================================================
% SECTION 5: Figure 2 - CAR vs power (mW)
% =========================================================================
fig2 = figure('Name', 'SPAD CAR vs Power', 'Color', 'w', 'Position', [80 100 800 480]);
hold on; box on; grid on; grid minor;

errorbar(powers_TE_mW, car_TE, err_car_TE, 'o-', ...
    'Color', c_TE, 'MarkerFaceColor', c_TE, 'MarkerSize', 7, ...
    'LineWidth', 1.5, 'DisplayName', 'TE');
errorbar(powers_TM_mW, car_TM, err_car_TM, 's-', ...
    'Color', c_TM, 'MarkerFaceColor', c_TM, 'MarkerSize', 7, ...
    'LineWidth', 1.5, 'DisplayName', 'TM');

set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickDir', 'out');
xlabel('Pump power (mW)', 'FontSize', 13);
ylabel('Coincidence-to-Accidental RatioAR', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 11);
hold off;
saveas(fig2, fullfile(outFolder, 'SPAD_CAR_vs_Power.png'));

% %% =========================================================================
% % SECTION 6: Figure 3 - CAR vs coincidence rate (true CC, cps)
% % =========================================================================
% % x-axis is the TRUE coincidence rate (the actual pair-generation rate),
% % since this is the physically meaningful "brightness vs quality"
% % trade-off curve. [PLACEHOLDER: swap to rawCC_TE/rawCC_TM below if you
% % instead want the axis to reflect total measured coincidence rate.]
% 
% fig3 = figure('Name', 'SPAD CAR vs Coincidence Rate', 'Color', 'w', 'Position', [80 100 800 480]);
% hold on; box on; grid on; grid minor;
% 
% errorbar(trueCC_TE, car_TE, err_car_TE, err_car_TE, err_trueCC_TE, err_trueCC_TE, 'o', ...
%     'Color', c_TE, 'MarkerFaceColor', c_TE, 'MarkerSize', 7, ...
%     'LineWidth', 1.5, 'DisplayName', 'TE');
% errorbar(trueCC_TM, car_TM, err_car_TM, err_car_TM, err_trueCC_TM, err_trueCC_TM, 's', ...
%     'Color', c_TM, 'MarkerFaceColor', c_TM, 'MarkerSize', 7, ...
%     'LineWidth', 1.5, 'DisplayName', 'TM');
% 
% set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickDir', 'out');
% xlabel('True coincidence rate (cps)', 'FontSize', 13);
% ylabel('CAR', 'FontSize', 13);
% legend('Location', 'northeast', 'FontSize', 11);
% hold off;
% saveas(fig3, fullfile(outFolder, 'SPAD_CAR_vs_CoincidenceRate.png'));

% %% =========================================================================
% % SECTION 7: Brightness - power-law fit to true coincidence rate vs power
% % =========================================================================
% % Brightness B is defined (Figures of Merit, Sec. 3.5.3) via
% %   R_true = B * P^2
% % in the low-power regime, before multipair generation degrades the
% % quadratic scaling. We identify "low-power" as all points up to and
% % including the CAR peak (beyond which accidentals/multipair effects
% % dominate and the simple quadratic model no longer applies).
% %
% % Two fits are reported:
% %   (1) FIXED exponent n=2, single-parameter fit through the origin:
% %       B_fixed = sum(P^2 .* R) / sum(P^4)   [matches thesis definition]
% %   (2) FREE exponent power-law fit (log-log linear regression), to
% %       check whether n really is close to 2.
% 
% [B_fixed_TE, n_free_TE, B_free_TE, nPts_TE] = ...
%     fitBrightness(powers_TE_mW, trueCC_TE, car_TE, 'TE');
% [B_fixed_TM, n_free_TM, B_free_TM, nPts_TM] = ...
%     fitBrightness(powers_TM_mW, trueCC_TM, car_TM, 'TM');
% 
% % --- Diagnostic log-log plot ---
% fig4 = figure('Name', 'Brightness fit diagnostic', 'Color', 'w', 'Position', [80 100 850 480]);
% hold on; box on; grid on;
% 
% loglog(powers_TE_mW, trueCC_TE, 'o', 'Color', c_TE_light, ...
%     'MarkerFaceColor', c_TE_light, 'MarkerSize', 6, 'DisplayName', 'TE (all points)');
% loglog(powers_TM_mW, trueCC_TM, 's', 'Color', c_TM_light, ...
%     'MarkerFaceColor', c_TM_light, 'MarkerSize', 6, 'DisplayName', 'TM (all points)');
% 
% loglog(powers_TE_mW(1:nPts_TE), trueCC_TE(1:nPts_TE), 'o', 'Color', c_TE, ...
%     'MarkerFaceColor', c_TE, 'MarkerSize', 7, 'DisplayName', 'TE (used in fit)');
% loglog(powers_TM_mW(1:nPts_TM), trueCC_TM(1:nPts_TM), 's', 'Color', c_TM, ...
%     'MarkerFaceColor', c_TM, 'MarkerSize', 7, 'DisplayName', 'TM (used in fit)');
% 
% xFitTE = powers_TE_mW(1:nPts_TE);
% xFitTM = powers_TM_mW(1:nPts_TM);
% loglog(xFitTE, B_free_TE * xFitTE.^n_free_TE, '-', 'Color', c_TE, ...
%     'LineWidth', 1.5, 'HandleVisibility', 'off');
% loglog(xFitTM, B_free_TM * xFitTM.^n_free_TM, '-', 'Color', c_TM, ...
%     'LineWidth', 1.5, 'HandleVisibility', 'off');
% 
% set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 11, 'TickDir', 'out');
% xlabel('Pump power (mW)', 'FontSize', 13);
% ylabel('True coincidence rate (cps)', 'FontSize', 13);
% legend('Location', 'northwest', 'FontSize', 9);
% hold off;
% saveas(fig4, fullfile(outFolder, 'Brightness_Fit_Diagnostic.png'));
% 
% fprintf('\nAll figures saved to: %s\n', outFolder);


%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function sigma_R = poissonRateError(R, T)
    R = max(R, 0);
    sigma_R = sqrt(R ./ T);
end


function cOut = lightenColor(cIn, factor)
% Blends a color toward white by `factor` (0 = unchanged, 1 = white).
    cOut = cIn + (1 - cIn) * factor;
end


function plotTransparentQuadFit(P, R, color)
% Overlays a quadratic fit (a*x^2 + b*x + c) in a transparent version of
% `color`, excluded from the legend.
    if length(P) < 4
        return
    end
    p = polyfit(P, R, 2);
    xFit = linspace(min(P), max(P), 200);
    yFit = polyval(p, xFit);

    h = plot(xFit, yFit, '-', 'Color', color, 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    try
        h.Color(4) = 0.35;   % add alpha transparency if supported
    catch
        % fallback: already using a light tint for accidentals; for the
        % darker true-CC fits, lighten further so it still reads as
        % "behind" the data points.
        h.Color = lightenColor(color, 0.5);
    end
end


function [B_fixed, n_free, B_free, nPts] = fitBrightness(P_mW, trueCC, CAR, label)
% Determines the low-power (pre-CAR-peak) subset, then fits brightness
% two ways: fixed exponent n=2 (through origin), and free exponent
% (log-log linear regression). Warns if too few points for a reliable
% fit.

    [~, peakIdx] = max(CAR);
    nPts = peakIdx;

    fprintf('\n--- Brightness fit: %s ---\n', label);
    fprintf('CAR peak at point %d of %d (P = %.2f mW)\n', peakIdx, length(P_mW), P_mW(peakIdx));

    if nPts < 4
        warning(['%s: only %d points before the CAR peak - brightness fit ' ...
                 'is likely unreliable. Consider reporting this qualitatively ' ...
                 'only, or using all available points with this caveat stated.'], ...
                 label, nPts);
    end

    P_low = P_mW(1:nPts);
    R_low = trueCC(1:nPts);

    % --- (1) Fixed exponent n=2, single-parameter, through origin ---
    B_fixed = sum((P_low.^2) .* R_low) / sum(P_low.^4);

    % --- (2) Free exponent, log-log linear regression ---
    % (guards against P=0 or R<=0 entries, which cannot be log-transformed)
    valid = (P_low > 0) & (R_low > 0);
    if sum(valid) >= 2
        logP = log10(P_low(valid));
        logR = log10(R_low(valid));
        p = polyfit(logP, logR, 1);
        n_free = p(1);
        B_free = 10^p(2);
    else
        n_free = NaN; B_free = NaN;
        warning('%s: not enough positive points for free-exponent fit.', label);
    end

    fprintf('Points used: %d\n', nPts);
    fprintf('Fixed-exponent (n=2) brightness: B = %.4g cps/mW^2\n', B_fixed);
    fprintf('Free-exponent fit: n = %.3f (theory expects ~2), B = %.4g cps/mW^n\n', ...
            n_free, B_free);
end