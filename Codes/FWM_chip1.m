% --------------------------------------------------------------
%% Power vs wvl (from OSA)
%% FWM Normalized Conversion Efficiency
% --------------------------------------------------------------

% ----------
%% SETTINGS
% ----------
clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%-------------
%%  LOAD DATA
%-------------
dataFolder = '../Tests/stimFWM/wvg1_TE1';
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);

sweepData = cell(numFiles, 1);
for k = 1:numFiles
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    sweepData{k} = readmatrix(fullFileName);
end

% ------------
%% SINGLE PLOT
% ------------
for i = [1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 20, 30, 35, 40, 50, 70, 72, 75, 80, 85, 90, 95, 100, 110, 120, 130, 135, 140, 146]
    sanityCheck = sweepData{i};
    fig = figure(i);
    plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');
    xlabel('Wavelength'); ylabel('Power (dBm)');
    title(sprintf('Sweep Data for File %d', i));
    grid on;
    i = i+1;
end

% ------------------------
%% PEAK TRACKING PIPELINE
% ------------------------
% Physical constants
c = 299792458;
wl_pump = 1.5496e-06;
pump_window = 0.5e-9;
search_window = 1.0e-9

% --- Thresholds ---
snr_threshold = 3;          % lowered from 20: real but weak peaks were being missed
OSA_floor     = -80;        % dBm: OSA fill value for out-of-range bins, must be excluded

% --- Debug: visualize noise window for these file indices ---
% Set to [] to skip. Add any file number you want to inspect.
debug_files = [];

% ------------------------------------------------------------------

genPower  = zeros(numFiles, 2);
seedPower = zeros(numFiles, 2);
pumpPower = zeros(numFiles, 1);

for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);

    if mean(x) > 1000
        x = x * 1e-9;
    end

    % --- Extract pump power ---
    is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
    if any(is_pump)
        pumpPower(k) = max(y(is_pump));
    else
        pumpPower(k) = NaN;
    end

    % --- Find seed position ---
    not_pump  = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
    x_noPump  = x(not_pump);
    y_noPump  = y(not_pump);

    [maxSeedPower, seedIdx] = max(y_noPump);
    wl_seed = x_noPump(seedIdx);
    seedPower(k, 1) = wl_seed;
    seedPower(k, 2) = maxSeedPower;

    % --- Reject seed = pump ---
    if abs(wl_seed - wl_pump) < 0.0015e-6
        genPower(k,1) = NaN;
        genPower(k,2) = NaN;
        continue;
    end

    % --- Calculate idler theoretical position ---
    w_pump       = c / wl_pump;
    w_seed       = c / wl_seed;
    w_idler      = 2*w_pump - w_seed;
    wl_idlerTheo = c / w_idler;

    %fprintf('File %02d | wl_idlerTheo=%.2f nm | OSA range: %.2f to %.2f nm\n', ...
    %k, wl_idlerTheo*1e9, x(1)*1e9, x(end)*1e9);

    % --- Search window around theoretical idler ---
    %int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window));
    int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window)) ...
    & (abs(x - wl_pump) > pump_window) ...   % exclude pump
    & (abs(x - wl_seed) > pump_window);      % exclude seed

    if any(int)
        [maxY, relIdx] = max(y(int));

        % --- Noise floor from flanking bands OUTSIDE the search window ---
        % FIX 1: use flanks, not the window itself (peak inflates median)
        % FIX 2: exclude pump wings (they inflate noise near pump)
        % FIX 3: exclude OSA fill value -210 dBm (out-of-range bins)
        lo_band = (x >= (wl_idlerTheo - 3*search_window)) & (x < (wl_idlerTheo - search_window)) ...
                & (abs(x - wl_pump) > pump_window) ...
                & (y > OSA_floor);
        hi_band = (x >  (wl_idlerTheo + search_window)) & (x <= (wl_idlerTheo + 3*search_window)) ...
                & (abs(x - wl_pump) > pump_window) ...
                & (y > OSA_floor);
        flank_y = y(lo_band | hi_band);

%
        if numel(flank_y) >= 5
            local_noise_floor = median(flank_y);
        else
            valid_int = y(int) > -200;
            if any(valid_int)
                local_noise_floor = median(y(int(valid_int)));
            else
                local_noise_floor = -80;  % hardcoded physical floor as last resort
            end
        end
        %

    %     if numel(flank_y) >= 5
    %         local_noise_floor = median(flank_y);
    %     else
    %         local_noise_floor = median(y(int));  % fallback if flanks are out of OSA range
    %     end
    %     fprintf('         maxY=%.1f dBm at %.2f nm | NF=%.1f dBm | SNR=%.1f dB\n', ...
    % maxY, x(find(int, 1))*1e9, local_noise_floor, maxY-local_noise_floor);
    %     % --- Debug plot for selected files ---
    %     if ismember(k, debug_files)
    %         figure('Name', sprintf('Debug File %d', k)); hold on;
    %         plot(x*1e9, y, 'Color', '#333333', 'LineWidth', 1);
    %         % Flank (noise) points in red
    %         plot(x(lo_band | hi_band)*1e9, y(lo_band | hi_band), ...
    %              'ro', 'MarkerFaceColor','r', 'MarkerSize', 5, 'DisplayName', 'Noise samples');
    %         % Search window points in blue
    %         plot(x(int)*1e9, y(int), ...
    %              'bs', 'MarkerFaceColor','b', 'MarkerSize', 5, 'DisplayName', 'Search window');
    %         xline(wl_idlerTheo*1e9, '--g', 'LineWidth', 1.5, 'Label', 'Idler theo', ...
    %               'LabelOrientation', 'horizontal');
    %         xline(wl_pump*1e9,      '--b', 'LineWidth', 1.5, 'Label', 'Pump', ...
    %               'LabelOrientation', 'horizontal');
    %         xline(wl_seed*1e9,      '--m', 'LineWidth', 1.5, 'Label', 'Seed', ...
    %               'LabelOrientation', 'horizontal');
    %         yline(local_noise_floor, '--r', 'LineWidth', 1.5, ...
    %               'Label', sprintf('NF = %.1f dBm', local_noise_floor), ...
    %               'LabelHorizontalAlignment', 'left');
    %         yline(maxY, '--k', 'LineWidth', 1, ...
    %               'Label', sprintf('Peak = %.1f dBm', maxY), ...
    %               'LabelHorizontalAlignment', 'left');
    %         xlabel('Wavelength (nm)'); ylabel('Power (dBm)'); grid on;
    %         title(sprintf('File %d $|$ SNR = %.1f dB $|$ threshold = %d dB', ...
    %                       k, maxY - local_noise_floor, snr_threshold));
    %         legend('Spectrum', 'Noise samples', 'Search window', ...
    %                'Location', 'best');
    %     end

        % --- SNR check ---
        if (maxY - local_noise_floor) >= snr_threshold
            idx = find(int);
            genPower(k, 1) = x(idx(relIdx));
            genPower(k, 2) = maxY;
        else
            genPower(k, 1) = NaN;
            genPower(k, 2) = NaN;
        end
    else
        genPower(k, 1) = NaN;
        genPower(k, 2) = NaN;
    end
end

% -----------------------------------------
%% COMPROVACIÓ MANUAL PER A L'OSA
% -----------------------------------------
disp('--- RESULTATS DEL PEAK TRACKING PER FITXER ---');
fprintf('SNR threshold: %d dB\n', snr_threshold);
for k = 1:numFiles
    if isnan(genPower(k,1))
        fprintf('Plot / Fitxer %02d: REBUTJAT (No supera el llindar SNR)\n', k);
    else
        fprintf('Plot / Fitxer %02d: Idler trobat a %.2f nm | Potència: %.2f dBm\n', ...
            k, genPower(k,1)*1e9, genPower(k,2));
    end
end
disp('----------------------------------------------');

% -------------------
%% NORMALIZED FWM VISUALIZATION
% -------------------
% 1. Clean up the data
valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);

wl_idler_valid = genPower(valid, 1) * 1e9;
P_idler_dBm    = genPower(valid, 2);
P_seed_dBm     = seedPower(valid, 2);
P_pump_dBm     = pumpPower(valid);

% Sort by wavelength to prevent zig-zag plots
[wl_idler_valid, sortIdx] = sort(wl_idler_valid);
P_idler_dBm = P_idler_dBm(sortIdx);
P_seed_dBm  = P_seed_dBm(sortIdx);
P_pump_dBm  = P_pump_dBm(sortIdx);

% 2. Convert dBm to Watts
P_idler_W = 10.^((P_idler_dBm - 30) / 10);
P_seed_W  = 10.^((P_seed_dBm  - 30) / 10);
P_pump_W  = 10.^((P_pump_dBm  - 30) / 10);

% 3. Normalized Efficiency
norm_eff_W2    = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));
norm_eff_dB_W2 = 10 * log10(norm_eff_W2);

% 4. Plot
fig = figure('Name', 'Normalized FWM Efficiency');
plot(wl_idler_valid, norm_eff_dB_W2, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
grid on;
xlabel('Generated photon wavelength (nm)');
ylabel('Normalized Efficiency $10\log_{10}(\eta_{norm})$ [dB W$^{-2}$]');
title('Normalized FWM Efficiency ($\eta_{norm}$)');
ylim([0 100])

% 5. Save
outFolder = '../Tests/stimFWM/wvg1_TE1';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig, fullfile(outFolder, 'Norm_ce.png'));


% -------------------
%% SEED POWER vs WAVELENGTH
% -------------------
fig_seed = figure('Name', 'Seed Power vs Wavelength');
plot(seedPower(:,1)*1e9, seedPower(:,2), '-squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
grid on;
xlabel('Seed wavelength (nm)');
ylabel('Seed Power (dBm)');
title('Seed Power vs Wavelength');
saveas(fig_seed, fullfile(outFolder, 'Seed_power.png'));