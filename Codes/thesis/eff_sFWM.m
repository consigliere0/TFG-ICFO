% =========================================================================
% TFG - ANALISI DE FOUR-WAVE MIXING (FWM) COMPARTIU (wvg3 vs wvg6)
% =========================================================================
% ÍNDEX DEL CODI:
%   SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
%   SECCIÓ 2: CONFIGURACIÓ DE RUTES I RECURSOS (Folders & Colors)
%   SECCIÓ 3: CANAL DE PROCESSAMENT AUTOMÀTIC (Pipeline Loop)
%   SECCIÓ 4: PLOT - COMPARATIVA D'EFICIÈNCIA NORMALITZADA (TE/TM Superposed)
% =========================================================================

% =========================================================================
%% SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
% =========================================================================
clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
% =========================================================================
%% SECCIÓ 2: CONFIGURACIÓ DE RUTES I RECURSOS (Folders & Colors)
% =========================================================================

% Només processem les dues polaritzacions de waveguide 2
dataFolders = { ...
    '../../Tests/stimFWM/wvg2_TE1', ... % Waveguide 2 TE
    '../../Tests/stimFWM/wvg2_TM1'  ... % Waveguide 2 TM
};

% Etiquetes per a la llegenda
setLabels = { ...
    'TE mode', ...
    'TM mode' ...
};

% Colors per a TE i TM de waveguide 2
setColors = { ...
    '#7E57C2'  ... % TM: lila fosc
    '#B39DDB', ... % TE: lila clar
    '#4CAF50'  ... % Loser TM: Verd pastel fosc
    '#A5D6A7', ... % Loser TE: Verd pastel clar
};

% Carpeta de sortida
outFolder = '../Tests/stimFWM/Comparison_Results';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

% Estructura on guardarem els resultats calculats
processedResults = cell(length(dataFolders), 1);

% =========================================================================
%% SECCIÓ 3: CANAL DE PROCESSAMENT AUTOMÀTIC (Pipeline Loop)
% =========================================================================
% --- Constants Físiques Generals ---
c = 299792458;
wl_pump = 1.5496e-06;
pump_window = 0.5e-9;
search_window = 1.0e-9;
snr_threshold = 3.2;        
OSA_floor     = -80;        

for f = 1:length(dataFolders)
    currentFolder = dataFolders{f};
    fprintf('Processant directori: %s...\n', currentFolder);

    filePattern = fullfile(currentFolder, '*.csv');
    csvFiles = dir(filePattern);
    numFiles = length(csvFiles);

    if numFiles == 0
        warning('No s''han trobat fitxers CSV a: %s. Saltant canal.', currentFolder);
        continue;
    end

    % Carrega de dades d'aquest directori
    sweepData = cell(numFiles, 1);
    for k = 1:numFiles
        baseFileName = csvFiles(k).name;
        fullFileName = fullfile(csvFiles(k).folder, baseFileName);
        sweepData{k} = readmatrix(fullFileName);
    end

    % Inicialització de matrius per al Peak Tracking d'aquest directori
    genPower  = zeros(numFiles, 2);
    seedPower = zeros(numFiles, 2);
    pumpPower = zeros(numFiles, 1);

    for k = 1:numFiles
        x = sweepData{k}(:,1);
        y = sweepData{k}(:,2);

        if mean(x) > 1000, x = x * 1e-9; end

        % Extreure Potència del Pump
        is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
        if any(is_pump)
            pumpPower(k) = max(y(is_pump));
        else
            pumpPower(k) = NaN;
        end

        % Extreure Potència i Posició del Seed
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

        % Càlcul de la posició teòrica de l'Idler
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

    % Neteja de dades i càlcul de normalització
    valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);
    wl_idler_valid = genPower(valid, 1) * 1e9;
    P_idler_dBm    = genPower(valid, 2);
    P_seed_dBm     = seedPower(valid, 2);
    P_pump_dBm     = pumpPower(valid);

    [wl_idler_valid, sortIdx] = sort(wl_idler_valid);
    P_idler_dBm = P_idler_dBm(sortIdx);
    P_seed_dBm  = P_seed_dBm(sortIdx);
    P_pump_dBm  = P_pump_dBm(sortIdx);

    P_idler_W = 10.^((P_idler_dBm - 30) / 10);
    P_seed_W  = 10.^((P_seed_dBm  - 30) / 10);
    P_pump_W  = 10.^((P_pump_dBm  - 30) / 10);

    norm_eff_W2    = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));
    norm_eff_dB_W2 = 10 * log10(norm_eff_W2);

    % Guardem les variables finals indexades per a fer el plot superposat després
    results.wl = wl_idler_valid;
    results.eff = norm_eff_dB_W2;
    processedResults{f} = results;
end

% =========================================================================
%% SECCIÓ 4: PLOT - EFICIÈNCIA NORMALITZADA DE WAVEGUIDE 2
% =========================================================================

fig_wvg2 = figure( ...
    'Name', 'FWM Normalized Efficiency - Waveguide 2', ...
    'Position', [100, 100, 850, 550]);

hold on;

plotsForLegend = [];

for f = 1:length(dataFolders)
    if isempty(processedResults{f})
        continue;
    end

    res = processedResults{f};

    % TE o TM de waveguide 2
    p = plot(res.wl, res.eff, '-s', ...
        'LineWidth', 1.8, ...
        'Color', setColors{f}, ...
        'MarkerFaceColor', setColors{f}, ...
        'MarkerSize', 5);

    plotsForLegend = [plotsForLegend, p];
end

grid on;

xlabel('Generated photon wavelength (nm)', 'FontSize', 11);
ylabel('Normalized Efficiency [dB W$^{-2}$]', ...
    'FontSize', 11);

%title('Sagnac Loop FWM Performance: $wvg_2$', 'FontSize', 12);

ylim([0 40]);

legend(plotsForLegend, setLabels, ...
    'Location', 'northeast', ...
    'FontSize', 10);

hold off;

% % Desem el gràfic només de waveguide 2
% outputFile = fullfile(outFolder, 'FWM_Efficiency_wvg2.png');
% saveas(fig_wvg2, outputFile);
% 
% fprintf('\nProcés completat! S''ha desat el gràfic de wvg2 a: %s\n', outputFile);
fprintf('\nFinished!\n')