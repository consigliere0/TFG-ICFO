% =========================================================================
% TFG - STIMULATED FWM: CENTRALIZED RESULTS PIPELINE
% =========================================================================
% ÍNDEX DEL CODI:
%   SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
%   SECCIÓ 2: CONFIGURACIÓ DE RUTES I RECURSOS (Folders & Colors)
%   SECCIÓ 3: BATCH PROCESSING (all 12 folders -> results struct array)
%   SECCIÓ 4: GAUSSIAN FITS (linear units) -> peak eff, peak wl, 3dB BW
%   SECCIÓ 5: FIGURE 1 - Representative single OSA spectrum (labeled)
%   SECCIÓ 6: FIGURE 2 - Spectral evolution overlay (winning device)
%   SECCIÓ 7: FIGURE 3 - Normalized efficiency, TE/TM, winner + Gaussian fit
%   SECCIÓ 8: FIGURE 4 - Bar chart, peak efficiency across all 6 devices
%   SECCIÓ 9: FIGURE 5 - Seed power ripple (Fabry-Perot artifact)
%   SECCIÓ 10: SUMMARY TABLE (six-device comparison) -> CSV export
%
% REQUIRES: Curve Fitting Toolbox (for `fit`, `fittype`, `confint`)
% =========================================================================

%% =========================================================================
% SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
% =========================================================================
clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% --- Physical constants and processing thresholds (from your original pipeline) ---
c              = 299792458;
wl_pump        = 1.5496e-06;   % [PLACEHOLDER: confirm exact pump wavelength used]
pump_window    = 0.5e-9;
search_window  = 1.0e-9;
snr_threshold  = 3.2;
OSA_floor      = -80;

%% =========================================================================
% SECCIÓ 2: CONFIGURACIÓ DE RUTES I RECURSOS (Folders & Colors)
% =========================================================================
% [PLACEHOLDER: confirm folder names/numbering match your actual six devices]
deviceNames = {'wvg1','wvg2','wvg3','wvg4','wvg5','wvg6'};
baseFolder  = '../../Tests/stimFWM/';

% Winning device used for Figures 1, 2, 3, 5 (single-device illustrative plots)
winnerDevice = 'wvg2';   % [PLACEHOLDER: confirm this is your selected device]

% Build the full 12-entry folder/label list (6 devices x TE/TM)
allFolders = {}; allLabels = {}; allDevice = {}; allPol = {};
for d = 1:length(deviceNames)
    allFolders{end+1} = fullfile(baseFolder, [deviceNames{d} '_TE1']); %#ok<*SAGROW>
    allLabels{end+1}  = sprintf('%s (TE)', deviceNames{d});
    allDevice{end+1}  = deviceNames{d};
    allPol{end+1}     = 'TE';

    allFolders{end+1} = fullfile(baseFolder, [deviceNames{d} '_TM1']);
    allLabels{end+1}  = sprintf('%s (TM)', deviceNames{d});
    allDevice{end+1}  = deviceNames{d};
    allPol{end+1}     = 'TM';
end

% Pastel color palette, extended from your original two-tone scheme.
% One light/dark pair per device (light = TE, dark = TM), 6 devices total.
devicePalette = { ...
    {'#B3D9FF', '#3A7CA5'}, ... % wvg1: light/dark blue
    {'#B39DDB', '#7E57C2'}, ... % wvg2: light/dark purple (your original "winner" pair)
    {'#A5D6A7', '#4CAF50'}, ... % wvg3: light/dark green   (your original "loser" pair)
    {'#FFCC80', '#EF6C00'}, ... % wvg4: light/dark orange
    {'#F48FB1', '#C2185B'}, ... % wvg5: light/dark pink
    {'#FFF59D', '#F9A825'}  ... % wvg6: light/dark yellow
};

outFolder = fullfile(baseFolder, 'Results_Figures');
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

%% =========================================================================
% SECCIÓ 3: BATCH PROCESSING (all 12 folders -> results struct array)
% =========================================================================
results = struct([]);

for f = 1:length(allFolders)
    fprintf('Processing %s ...\n', allFolders{f});
    R = processFolder(allFolders{f}, c, wl_pump, pump_window, ...
                       search_window, snr_threshold, OSA_floor);
    R.label  = allLabels{f};
    R.device = allDevice{f};
    R.pol    = allPol{f};
    R.folder = allFolders{f};
    results = [results, R]; %#ok<AGROW>
end

%% =========================================================================
% SECCIÓ 4: GAUSSIAN FITS (linear units) -> peak eff, peak wl, 3dB BW
% =========================================================================
% IMPORTANT: fit is performed on LINEAR efficiency (norm_eff_W2), never on
% the dB-scaled quantity, since FWHM is only physically meaningful in
% linear power units (see thesis discussion).

for f = 1:length(results)
    wl  = results(f).wl_idler_valid;   % [nm]
    eff = results(f).norm_eff_W2;       % linear, [W^-2]

    if numel(wl) < 4
        warning('Too few points to fit %s - skipping.', results(f).label);
        results(f).fit = [];
        results(f).peakEff_dB = NaN;
        results(f).peakWl = NaN;
        results(f).bw3dB_nm = NaN;
        results(f).bw3dB_CI = [NaN NaN];
        continue
    end

    ft = fittype('gauss1');
    try
        [fr, gof] = fit(wl, eff, ft); %#ok<ASGLU>
        FWHM_nm = 2 * fr.c1 * sqrt(log(2));
        ci = confint(fr, 0.95);
        bw_CI = 2*sqrt(log(2)) * [ci(1,3), ci(2,3)];

        results(f).fit         = fr;
        results(f).peakEff_dB  = 10*log10(fr.a1);
        results(f).peakWl      = fr.b1;
        results(f).bw3dB_nm    = FWHM_nm;
        results(f).bw3dB_CI    = bw_CI;
    catch ME
        warning('Fit failed for %s: %s', results(f).label, ME.message);
        results(f).fit = [];
        results(f).peakEff_dB = NaN;
        results(f).peakWl = NaN;
        results(f).bw3dB_nm = NaN;
        results(f).bw3dB_CI = [NaN NaN];
    end
end

%% =========================================================================
% SECCIÓ 5: FIGURE 1 - Representative single OSA spectrum (labeled)
% =========================================================================
% Pick one clean, high-SNR file from the winning device (TE) to illustrate
% a single raw measurement, with pump/seed/idler peaks labeled.

winIdx_TE = find(strcmp({results.device}, winnerDevice) & strcmp({results.pol}, 'TE'), 1);
Rwin = results(winIdx_TE);

exampleFileIdx = 100;   % [PLACEHOLDER: pick a specific file index with clean, high-SNR peaks]
if exampleFileIdx <= numel(Rwin.sweepData)
    x = Rwin.sweepData{exampleFileIdx}(:,1);
    y = Rwin.sweepData{exampleFileIdx}(:,2);
    if mean(x) < 1000, x = x * 1e9; end  % ensure nm

    fig1 = figure('Name', 'Representative OSA Spectrum', 'Position', [100 100 850 500]);
    plot(x, y, 'Color', '#546E7A', 'LineWidth', 1.3);
    grid on; box on;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 11);
    xlabel('Wavelength (nm)', 'FontSize', 12);
    ylabel('Power (dBm)', 'FontSize', 12);
    title(sprintf('Representative stimulated FWM spectrum (%s)', Rwin.label), 'FontSize', 13);

    % Annotate pump, seed, idler positions if within this trace's valid data
    [~, pumpIdx] = min(abs(x - wl_pump*1e9));
    text(x(pumpIdx), y(pumpIdx)+3, 'Pump', 'FontSize', 10, 'HorizontalAlignment', 'center');

    fIdxInValid = find(Rwin.fileIdx_valid == exampleFileIdx, 1);
    if ~isempty(fIdxInValid)
        seed_wl = Rwin.seedPower_all(exampleFileIdx,1)*1e9;
        seed_pw = Rwin.seedPower_all(exampleFileIdx,2);
        idler_wl = Rwin.wl_idler_valid(fIdxInValid);
        idler_pw = Rwin.P_idler_dBm(fIdxInValid);
        text(seed_wl, seed_pw+3, 'Seed', 'FontSize', 10, 'HorizontalAlignment', 'center');
        text(idler_wl, idler_pw+3, 'Generated idler', 'FontSize', 10, 'HorizontalAlignment', 'center');
    end

    saveas(fig1, fullfile(outFolder, 'Fig1_Representative_Spectrum.png'));
else
    warning('exampleFileIdx out of range - skipping Figure 1.');
end

% %% =========================================================================
% % SECCIÓ 6: FIGURE 2 - Spectral evolution overlay (winning device)
% % =========================================================================
% fig2 = figure('Name', 'Spectra Overlay Publication', 'Position', [150, 150, 950, 520]);
% hold on;
% 
% subset_files = [1, 5, 10, 30, 80, 100, 120, 140];   % [PLACEHOLDER: adjust to your file count]
% nLines = length(subset_files);
% 
% c1 = [0.12, 0.47, 0.71];  % deep blue
% c2 = [0.73, 0.15, 0.48];  % magenta
% custom_map = [linspace(c1(1),c2(1),nLines)', linspace(c1(2),c2(2),nLines)', linspace(c1(3),c2(3),nLines)'];
% 
% for j = 1:nLines
%     idx = subset_files(j);
%     if idx <= numel(Rwin.sweepData)
%         x_nm = Rwin.sweepData{idx}(:,1);
%         y_dBm = Rwin.sweepData{idx}(:,2);
%         if mean(x_nm) < 1000, x_nm = x_nm * 1e9; end
% 
%         wl_seed_actual_nm = Rwin.seedPower_all(idx,1) * 1e9;
%         if ~isnan(wl_seed_actual_nm)
%             seed_mask = (x_nm >= (wl_seed_actual_nm - 1.5)) & (x_nm <= (wl_seed_actual_nm + 1.5));
%             y_dBm(seed_mask) = NaN;
%         end
% 
%         plot(x_nm, y_dBm, 'Color', custom_map(j,:), 'LineWidth', 1.3);
%     end
% end
% 
% xline(wl_pump*1e9, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8, ...
%       'Label', 'Pump', 'LabelOrientation', 'horizontal', 'FontSize', 10);
% 
% xlabel('Wavelength (nm)', 'FontSize', 12);
% ylabel('Power (dBm)', 'FontSize', 12);
% title(sprintf('Spectral evolution of the generated signal (%s)', Rwin.label), ...
%       'FontSize', 13, 'FontWeight', 'bold');
% 
% box on;
% set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11, ...
%     'GridColor', [0.85 0.85 0.85], 'GridAlpha', 0.5);
% grid on;
% ylim([-90 -20]);
% 
% colormap(custom_map);
% cb = colorbar;
% cb.Label.Interpreter = 'latex';
% cb.Label.String = 'Experimental sweep index ($\rightarrow$ seed wavelength)';
% cb.Label.FontSize = 11;
% set(cb, 'Ticks', [0, 0.5, 1], ...
%     'TickLabels', {num2str(subset_files(1)), num2str(median(subset_files)), num2str(subset_files(end))});
% 
% hold off;
% saveas(fig2, fullfile(outFolder, 'Fig2_Spectral_Overlay.png'));

%% =========================================================================
% SECCIÓ 7: FIGURE 3 - Normalized efficiency, TE/TM, winner + Gaussian fit
% =========================================================================
winIdx_TM = find(strcmp({results.device}, winnerDevice) & strcmp({results.pol}, 'TM'), 1);
Rte = results(winIdx_TE);
Rtm = results(winIdx_TM);

fig3 = figure('Name', 'FWM Normalized Efficiency - Winner TE/TM', 'Position', [100, 100, 850, 550]);
hold on;

% --- Raw data points ---
pTE = plot(Rte.wl_idler_valid, Rte.norm_eff_dB_W2, 'o', ...
    'MarkerSize', 5, 'MarkerFaceColor', '#B39DDB', 'MarkerEdgeColor', '#7E57C2');
pTM = plot(Rtm.wl_idler_valid, Rtm.norm_eff_dB_W2, 's', ...
    'MarkerSize', 5, 'MarkerFaceColor', '#A5D6A7', 'MarkerEdgeColor', '#4CAF50');

% --- Gaussian fit overlays (converted back to dB for display) + FWHM annotation ---
plotFitWithBandwidth(Rte, '#7E57C2');
plotFitWithBandwidth(Rtm, '#4CAF50');

grid on; box on;
set(gca, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 11);
xlabel('Generated photon wavelength (nm)', 'FontSize', 11);
ylabel('Normalized efficiency $10\log_{10}(\eta_{norm})$ [dB W$^{-2}$]', 'FontSize', 11);
title(sprintf('Stimulated FWM efficiency: %s', winnerDevice), 'FontSize', 12);
legend([pTE, pTM], {'TE data', 'TM data'}, 'Location', 'best', 'FontSize', 10);
ylim([0 100]);

hold off;
saveas(fig3, fullfile(outFolder, 'Fig3_Efficiency_Winner_TE_TM.png'));

fprintf('\n--- Winning device (%s) fit results ---\n', winnerDevice);
fprintf('TE: peak %.2f dB W^-2 at %.2f nm, 3-dB BW = %.2f nm [%.2f, %.2f]\n', ...
    Rte.peakEff_dB, Rte.peakWl, Rte.bw3dB_nm, Rte.bw3dB_CI(1), Rte.bw3dB_CI(2));
fprintf('TM: peak %.2f dB W^-2 at %.2f nm, 3-dB BW = %.2f nm [%.2f, %.2f]\n', ...
    Rtm.peakEff_dB, Rtm.peakWl, Rtm.bw3dB_nm, Rtm.bw3dB_CI(1), Rtm.bw3dB_CI(2));

% %% =========================================================================
% % SECCIÓ 8: FIGURE 4 - Bar chart, peak efficiency across all 6 devices
% % =========================================================================
% peakEff_TE = nan(1,6); peakEff_TM = nan(1,6);
% for d = 1:6
%     idxTE = find(strcmp({results.device}, deviceNames{d}) & strcmp({results.pol}, 'TE'), 1);
%     idxTM = find(strcmp({results.device}, deviceNames{d}) & strcmp({results.pol}, 'TM'), 1);
%     if ~isempty(idxTE), peakEff_TE(d) = results(idxTE).peakEff_dB; end
%     if ~isempty(idxTM), peakEff_TM(d) = results(idxTM).peakEff_dB; end
% end
% 
% fig4 = figure('Name', 'Peak Efficiency Comparison', 'Position', [100 100 850 500]);
% b = bar([peakEff_TE; peakEff_TM]');
% b(1).FaceColor = '#B39DDB'; b(1).EdgeColor = '#7E57C2'; b(1).LineWidth = 1.1;
% b(2).FaceColor = '#A5D6A7'; b(2).EdgeColor = '#4CAF50'; b(2).LineWidth = 1.1;
% 
% set(gca, 'XTickLabel', deviceNames, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 11);
% xlabel('Sagnac loop device', 'FontSize', 12);
% ylabel('Peak normalized efficiency [dB W$^{-2}$]', 'FontSize', 12);
% title('Peak stimulated FWM efficiency across all characterized devices', 'FontSize', 13);
% legend({'TE', 'TM'}, 'Location', 'best', 'FontSize', 10);
% grid on; box on;
% 
% saveas(fig4, fullfile(outFolder, 'Fig4_PeakEfficiency_AllDevices.png'));

%% =========================================================================
% SECCIÓ 9: FIGURE 5 - Seed power ripple (Fabry-Perot artifact)
% =========================================================================
fig5 = figure('Name', 'Seed Power Ripple', 'Position', [100 100 800 400]);
plot(Rte.clean_seed_wl_nm, Rte.clean_seed_pwr, '-s', 'LineWidth', 1.5, ...
     'Color', '#90CAF9', 'MarkerFaceColor', '#90CAF9', 'MarkerSize', 4);
grid on; box on;
set(gca, 'TickDir', 'out', 'LineWidth', 1.1, 'FontSize', 11);
xlabel('Seed wavelength (nm)', 'FontSize', 11);
ylabel('Seed power (dBm)', 'FontSize', 11);
title(sprintf('Input seed power sweep, showing Fabry-P\\''erot ripple (%s)', Rte.label), 'FontSize', 12);

saveas(fig5, fullfile(outFolder, 'Fig5_Seed_Power_Ripple.png'));

% %% =========================================================================
% % SECCIÓ 10: SUMMARY TABLE (six-device comparison) -> CSV export
% % =========================================================================
% Device = {results.device}';
% Polarization = {results.pol}';
% PeakEfficiency_dBW2 = [results.peakEff_dB]';
% PeakWavelength_nm = [results.peakWl]';
% Bandwidth3dB_nm = [results.bw3dB_nm]';
% BW_CI_low = cellfun(@(c) c(1), {results.bw3dB_CI})';
% BW_CI_high = cellfun(@(c) c(2), {results.bw3dB_CI})';
% 
% T = table(Device, Polarization, PeakEfficiency_dBW2, PeakWavelength_nm, ...
%           Bandwidth3dB_nm, BW_CI_low, BW_CI_high);
% 
% disp('=== Six-device stimulated FWM summary ===');
% disp(T);
% 
% writetable(T, fullfile(outFolder, 'stimFWM_summary_table.csv'));
% fprintf('\nSummary table written to: %s\n', fullfile(outFolder, 'stimFWM_summary_table.csv'));
% fprintf('All figures saved to: %s\n', outFolder);


%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function R = processFolder(dataFolder, c, wl_pump, pump_window, search_window, snr_threshold, OSA_floor)
% Loads all CSVs in dataFolder and runs the peak-tracking pipeline
% (pump / seed / generated-idler extraction + normalized efficiency),
% exactly as in the original per-folder script, but packaged for reuse.

    filePattern = fullfile(dataFolder, '*.csv');
    csvFiles = dir(filePattern);
    numFiles = length(csvFiles);

    R = struct();
    R.sweepData = cell(numFiles, 1);

    if numFiles == 0
        warning('No CSV files found in: %s', dataFolder);
        R.wl_idler_valid = []; R.norm_eff_W2 = []; R.norm_eff_dB_W2 = [];
        R.P_idler_dBm = []; R.seedPower_all = []; R.fileIdx_valid = [];
        R.clean_seed_wl_nm = []; R.clean_seed_pwr = [];
        return
    end

    for k = 1:numFiles
        fullFileName = fullfile(csvFiles(k).folder, csvFiles(k).name);
        R.sweepData{k} = readmatrix(fullFileName);
    end

    genPower  = zeros(numFiles, 2);
    seedPower = zeros(numFiles, 2);
    pumpPower = zeros(numFiles, 1);

    for k = 1:numFiles
        x = R.sweepData{k}(:,1);
        y = R.sweepData{k}(:,2);
        if mean(x) > 1000, x = x * 1e-9; end

        % Pump
        is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
        if any(is_pump)
            pumpPower(k) = max(y(is_pump));
        else
            pumpPower(k) = NaN;
        end

        % Seed
        not_pump  = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
        x_noPump  = x(not_pump);
        y_noPump  = y(not_pump);
        [maxSeedPower, seedIdx] = max(y_noPump);
        wl_seed = x_noPump(seedIdx);
        seedPower(k, 1) = wl_seed;
        seedPower(k, 2) = maxSeedPower;

        if abs(wl_seed - wl_pump) < 0.0015e-6
            genPower(k,1) = NaN; genPower(k,2) = NaN;
            continue;
        end

        % Theoretical idler position
        w_pump       = c / wl_pump;
        w_seed       = c / wl_seed;
        w_idler      = 2*w_pump - w_seed;
        wl_idlerTheo = c / w_idler;

        int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window)) ...
            & (abs(x - wl_pump) > pump_window) ...
            & (abs(x - wl_seed) > pump_window);

        if any(int)
            [maxY, relIdx] = max(y(int));

            lo_band = (x >= (wl_idlerTheo - 3*search_window)) & (x < (wl_idlerTheo - search_window)) ...
                    & (abs(x - wl_pump) > pump_window) & (y > OSA_floor);
            hi_band = (x >  (wl_idlerTheo + search_window)) & (x <= (wl_idlerTheo + 3*search_window)) ...
                    & (abs(x - wl_pump) > pump_window) & (y > OSA_floor);
            flank_y = y(lo_band | hi_band);

            if numel(flank_y) >= 5
                local_noise_floor = median(flank_y);
            else
                valid_int = y(int) > -200;
                if any(valid_int)
                    local_noise_floor = median(y(int(valid_int)));
                else
                    local_noise_floor = -80;
                end
            end

            if (maxY - local_noise_floor) >= snr_threshold
                idx = find(int);
                genPower(k, 1) = x(idx(relIdx));
                genPower(k, 2) = maxY;
            else
                genPower(k, 1) = NaN; genPower(k, 2) = NaN;
            end
        else
            genPower(k, 1) = NaN; genPower(k, 2) = NaN;
        end
    end

    valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);
    fileIdx_all = (1:numFiles)';

    wl_idler_valid = genPower(valid, 1) * 1e9;
    P_idler_dBm    = genPower(valid, 2);
    P_seed_dBm     = seedPower(valid, 2);
    P_pump_dBm     = pumpPower(valid);
    fileIdx_valid  = fileIdx_all(valid);

    [wl_idler_valid, sortIdx] = sort(wl_idler_valid);
    P_idler_dBm   = P_idler_dBm(sortIdx);
    P_seed_dBm    = P_seed_dBm(sortIdx);
    P_pump_dBm    = P_pump_dBm(sortIdx);
    fileIdx_valid = fileIdx_valid(sortIdx);

    P_idler_W = 10.^((P_idler_dBm - 30) / 10);
    P_seed_W  = 10.^((P_seed_dBm  - 30) / 10);
    P_pump_W  = 10.^((P_pump_dBm  - 30) / 10);

    norm_eff_W2    = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));
    norm_eff_dB_W2 = 10 * log10(norm_eff_W2);

    % Seed power trace with pump-overlap point removed (for ripple plot)
    seed_wl_nm = seedPower(:,1) * 1e9;
    seed_pwr   = seedPower(:,2);
    valid_seed_idx = abs(seed_wl_nm - wl_pump*1e9) > 0.5 & ~isnan(seed_wl_nm);

    R.wl_idler_valid   = wl_idler_valid;
    R.P_idler_dBm      = P_idler_dBm;
    R.norm_eff_W2      = norm_eff_W2;
    R.norm_eff_dB_W2   = norm_eff_dB_W2;
    R.seedPower_all    = seedPower;      % [wl(m), power(dBm)] per file, unfiltered
    R.fileIdx_valid    = fileIdx_valid;
    R.clean_seed_wl_nm = seed_wl_nm(valid_seed_idx);
    R.clean_seed_pwr   = seed_pwr(valid_seed_idx);
end


function plotFitWithBandwidth(R, colorHex)
% Overlays a Gaussian fit (converted to dB) on the current axes, and
% annotates the 3-dB bandwidth with a horizontal double-headed arrow.

    if isempty(R.fit) || isnan(R.bw3dB_nm)
        return
    end

    xFit = linspace(min(R.wl_idler_valid), max(R.wl_idler_valid), 500)';
    yFit_linear = R.fit(xFit);
    yFit_dB = 10*log10(yFit_linear);

    plot(xFit, yFit_dB, '-', 'Color', colorHex, 'LineWidth', 1.8);

    % Half-max (3 dB down) level and FWHM span
    halfMax_dB = R.peakEff_dB - 3.01;
    wl_lo = R.peakWl - R.bw3dB_nm/2;
    wl_hi = R.peakWl + R.bw3dB_nm/2;

    plot([wl_lo wl_hi], [halfMax_dB halfMax_dB], '--', 'Color', colorHex, 'LineWidth', 1);
    plot([wl_lo wl_lo], [halfMax_dB-3, halfMax_dB+3], ':', 'Color', colorHex, 'LineWidth', 0.8);
    plot([wl_hi wl_hi], [halfMax_dB-3, halfMax_dB+3], ':', 'Color', colorHex, 'LineWidth', 0.8);

    text(R.peakWl, halfMax_dB - 6, sprintf('BW$_{3dB}$ = %.1f nm', R.bw3dB_nm), ...
        'Color', colorHex, 'FontSize', 9, 'HorizontalAlignment', 'center');
end