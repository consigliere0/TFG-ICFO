% % =========================================================================
% % TFG - ANALISI DE FOUR-WAVE MIXING (FWM) COMPARTIU (wvg3 vs wvg6)
% % =========================================================================
% % ÍNDEX DEL CODI:
% %   SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
% %   SECCIÓ 2: CONFIGURACIÓ DE RUTES I RECURSOS (Folders & Colors)
% %   SECCIÓ 3: CANAL DE PROCESSAMENT AUTOMÀTIC (Pipeline Loop)
% %   SECCIÓ 4: PLOT - COMPARATIVA D'EFICIÈNCIA NORMALITZADA (TE/TM Superposed)
% % =========================================================================
% 
% % =========================================================================
% %% SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
% % =========================================================================
% clear; close all; format long g
% s = settings;
% s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
% set(groot, 'defaultTextInterpreter', 'latex')
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
% set(groot, 'defaultLegendInterpreter','latex');
% 
% % =========================================================================
% %% SECCIÓ 2: CONFIGURACIÓ DE RUTES I RECURSOS (Folders & Colors)
% % =========================================================================
% % Definim els directoris dels conjunts que volem comparar
% dataFolders = { ...
%     '../Tests/stimFWM/wvg2_TE1', ... % Guanyador TE
%     '../Tests/stimFWM/wvg2_TM1', ... % Guanyador TM
%     '../Tests/stimFWM/wvg6_TE1', ... % Perdedor TE
%     '../Tests/stimFWM/wvg6_TM1'  ... % Perdedor TM
% };
% 
% % Etiquetes per a la llegenda del gràfic
% setLabels = { ...
%     'Winning $wvg_2$ (TE)', ...
%     'Winning $wvg_2$ (TM)', ...
%     'Loser $wvg_6$ (TE)', ...
%     'Loser $wvg_6$ (TM)' ...
% };
% 
% % Colors Pastel demanats (Dos tons de lila/purple, dos tons de verd/green)
% setColors = { ...
%     '#B39DDB', ... % Winning TE: Lila pastel clar
%     '#7E57C2', ... % Winning TM: Lila pastel fosc
%     '#A5D6A7', ... % Loser TE: Verd pastel clar
%     '#4CAF50'  ... % Loser TM: Verd pastel fosc
% };
% 
% % Carpeta de sortida dels gràfics comparatius
% outFolder = '../Tests/stimFWM/Comparison_Results';
% if ~exist(outFolder, 'dir'), mkdir(outFolder); end
% 
% % Estructura on guardarem els resultats calculats de cada canal
% processedResults = cell(length(dataFolders), 1);
% 
% % =========================================================================
% %% SECCIÓ 3: CANAL DE PROCESSAMENT AUTOMÀTIC (Pipeline Loop)
% % =========================================================================
% % --- Constants Físiques Generals ---
% c = 299792458;
% wl_pump = 1.5496e-06;
% pump_window = 0.5e-9;
% search_window = 1.0e-9;
% snr_threshold = 3.2;        
% OSA_floor     = -80;        
% 
% for f = 1:length(dataFolders)
%     currentFolder = dataFolders{f};
%     fprintf('Processant directori: %s...\n', currentFolder);
% 
%     filePattern = fullfile(currentFolder, '*.csv');
%     csvFiles = dir(filePattern);
%     numFiles = length(csvFiles);
% 
%     if numFiles == 0
%         warning('No s''han trobat fitxers CSV a: %s. Saltant canal.', currentFolder);
%         continue;
%     end
% 
%     % Carrega de dades d'aquest directori
%     sweepData = cell(numFiles, 1);
%     for k = 1:numFiles
%         baseFileName = csvFiles(k).name;
%         fullFileName = fullfile(csvFiles(k).folder, baseFileName);
%         sweepData{k} = readmatrix(fullFileName);
%     end
% 
%     % Inicialització de matrius per al Peak Tracking d'aquest directori
%     genPower  = zeros(numFiles, 2);
%     seedPower = zeros(numFiles, 2);
%     pumpPower = zeros(numFiles, 1);
% 
%     for k = 1:numFiles
%         x = sweepData{k}(:,1);
%         y = sweepData{k}(:,2);
% 
%         if mean(x) > 1000, x = x * 1e-9; end
% 
%         % Extreure Potència del Pump
%         is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
%         if any(is_pump)
%             pumpPower(k) = max(y(is_pump));
%         else
%             pumpPower(k) = NaN;
%         end
% 
%         % Extreure Potència i Posició del Seed
%         not_pump  = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
%         x_noPump  = x(not_pump);
%         y_noPump  = y(not_pump);
%         [maxSeedPower, seedIdx] = max(y_noPump);
%         wl_seed = x_noPump(seedIdx);
%         seedPower(k, 1) = wl_seed;
%         seedPower(k, 2) = maxSeedPower;
% 
%         if abs(wl_seed - wl_pump) < 0.0015e-6
%             genPower(k,1) = NaN; genPower(k,2) = NaN;
%             continue;
%         end
% 
%         % Càlcul de la posició teòrica de l'Idler
%         w_pump       = c / wl_pump;
%         w_seed       = c / wl_seed;
%         w_idler      = 2*w_pump - w_seed;
%         wl_idlerTheo = c / w_idler;
% 
%         int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window)) ...
%             & (abs(x - wl_pump) > pump_window) ...   
%             & (abs(x - wl_seed) > pump_window);      
% 
%         if any(int)
%             [maxY, relIdx] = max(y(int));
% 
%             lo_band = (x >= (wl_idlerTheo - 3*search_window)) & (x < (wl_idlerTheo - search_window)) ...
%                     & (abs(x - wl_pump) > pump_window) & (y > OSA_floor);
%             hi_band = (x >  (wl_idlerTheo + search_window)) & (x <= (wl_idlerTheo + 3*search_window)) ...
%                     & (abs(x - wl_pump) > pump_window) & (y > OSA_floor);
%             flank_y = y(lo_band | hi_band);
% 
%             if numel(flank_y) >= 5
%                 local_noise_floor = median(flank_y);
%             else
%                 valid_int = y(int) > -200;
%                 if any(valid_int)
%                     local_noise_floor = median(y(int(valid_int)));
%                 else
%                     local_noise_floor = -80; 
%                 end
%             end
% 
%             if (maxY - local_noise_floor) >= snr_threshold
%                 idx = find(int);
%                 genPower(k, 1) = x(idx(relIdx));
%                 genPower(k, 2) = maxY;
%             else
%                 genPower(k, 1) = NaN; genPower(k, 2) = NaN;
%             end
%         else
%             genPower(k, 1) = NaN; genPower(k, 2) = NaN;
%         end
%     end
% 
%     % Neteja de dades i càlcul de normalització
%     valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);
%     wl_idler_valid = genPower(valid, 1) * 1e9;
%     P_idler_dBm    = genPower(valid, 2);
%     P_seed_dBm     = seedPower(valid, 2);
%     P_pump_dBm     = pumpPower(valid);
% 
%     [wl_idler_valid, sortIdx] = sort(wl_idler_valid);
%     P_idler_dBm = P_idler_dBm(sortIdx);
%     P_seed_dBm  = P_seed_dBm(sortIdx);
%     P_pump_dBm  = P_pump_dBm(sortIdx);
% 
%     P_idler_W = 10.^((P_idler_dBm - 30) / 10);
%     P_seed_W  = 10.^((P_seed_dBm  - 30) / 10);
%     P_pump_W  = 10.^((P_pump_dBm  - 30) / 10);
% 
%     norm_eff_W2    = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));
%     norm_eff_dB_W2 = 10 * log10(norm_eff_W2);
% 
%     % Guardem les variables finals indexades per a fer el plot superposat després
%     results.wl = wl_idler_valid;
%     results.eff = norm_eff_dB_W2;
%     processedResults{f} = results;
% end
% 
% % =========================================================================
% %% SECCIÓ 4: PLOT - COMPARATIVA D'EFICIÈNCIA NORMALITZADA (Superposed Plots)
% % =========================================================================
% fig_comparison = figure('Name', 'FWM Normalized Efficiency Comparison', 'Position', [100, 100, 850, 550]);
% hold on;
% 
% plotsForLegend = []; % Vector temporal per gestionar la llegenda de colors netament
% 
% for f = 1:length(dataFolders)
%     if isempty(processedResults{f}), continue; end
% 
%     res = processedResults{f};
% 
%     % Dibuixem cada línia amb el seu respectiu color pastel de la configuració
%     % S'utilitzen quadrats per al xip ('square') per mantenir l'estil formal
%     p = plot(res.wl, res.eff, '-s', 'LineWidth', 1.8, ...
%              'Color', setColors{f}, ...
%              'MarkerFaceColor', setColors{f}, ...
%              'MarkerSize', 5);
% 
%     plotsForLegend = [plotsForLegend, p]; % Assignem l'element per a la llegenda
% end
% 
% grid on;
% xlabel('Generated photon wavelength (nm)', 'FontSize', 11);
% ylabel('Normalized Efficiency $10\log_{10}(\eta_{norm})$ [dB W$^{-2}$]', 'FontSize', 11);
% title('Sagnac Loop FWM Performance: $wvg_2$ (Winner) vs $wvg_6$ (Loser)', 'FontSize', 12);
% 
% % Límits eix Y configurats a [0 100] tal com tenies predefinit, modifica si canvia el rang
% ylim([0 100]); 
% 
% % Afegim la llegenda utilitzant intèrpret de LaTeX per admetre subíndexs ($wvg_3$)
% legend(plotsForLegend, setLabels, 'Location', 'best', 'FontSize', 10);
% 
% hold off;
% 
% % Desem el gràfic final a la carpeta compartida de resultats
% saveas(fig_comparison, fullfile(outFolder, 'FWM_Efficiency_Comparison.png'));
% fprintf('\nProcés completat! S''ha desat el gràfic comparatiu a: %s\n', fullfile(outFolder, 'FWM_Efficiency_Comparison.png'));



























%%
% % =========================================================================
% % TFG - AMITJÀ DE TESI: SUBPLOTS DE POTÈNCIA (Waveguide 2 TE)
% % =========================================================================
% % ÍNDEX DEL CODI:
% %   SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
% %   SECCIÓ 2: CÀRREGA DE DADES (Load Data)
% %   SECCIÓ 3: EXTRACCIÓ DE PICS (Peak Tracking Pipeline)
% %   SECCIÓ 4: NETEJA I ORDENACIÓ DE DADES
% %   SECCIÓ 5: FIGURA PÒSTER - SUBPLOTS SEED I GENERATED SIGNAL (Pastel)
% % =========================================================================
% 
% % =========================================================================
% %% SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
% % =========================================================================
% clear; close all; format long g
% s = settings;
% s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
% set(groot, 'defaultTextInterpreter', 'latex')
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
% set(groot, 'defaultLegendInterpreter','latex');
% 
% % =========================================================================
% %% SECCIÓ 2: CÀRREGA DE DADES (Load Data)
% % =========================================================================
% dataFolder = '../Tests/stimFWM/wvg2_TE1';
% filePattern = fullfile(dataFolder, '*.csv');
% csvFiles = dir(filePattern);
% numFiles = length(csvFiles);
% sweepData = cell(numFiles, 1);
% 
% for k = 1:numFiles
%     baseFileName = csvFiles(k).name;
%     fullFileName = fullfile(csvFiles(k).folder, baseFileName);
%     sweepData{k} = readmatrix(fullFileName);
% end
% 
% outFolder = dataFolder; 
% if ~exist(outFolder, 'dir'), mkdir(outFolder); end
% 
% % =========================================================================
% %% SECCIÓ 3: EXTRACCIÓ DE PICS (Peak Tracking Pipeline)
% % =========================================================================
% c = 299792458;
% wl_pump = 1.5496e-06;
% pump_window = 0.5e-9;
% search_window = 1.0e-9;
% snr_threshold = 3.2;        
% OSA_floor     = -80;        
% 
% genPower  = zeros(numFiles, 2);
% seedPower = zeros(numFiles, 2);
% pumpPower = zeros(numFiles, 1);
% 
% for k = 1:numFiles
%     x = sweepData{k}(:,1);
%     y = sweepData{k}(:,2);
% 
%     if mean(x) > 1000, x = x * 1e-9; end
% 
%     % Extreure Pump
%     is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
%     if any(is_pump)
%         pumpPower(k) = max(y(is_pump));
%     else
%         pumpPower(k) = NaN;
%     end
% 
%     % Extreure Seed
%     not_pump  = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
%     x_noPump  = x(not_pump);
%     y_noPump  = y(not_pump);
%     [maxSeedPower, seedIdx] = max(y_noPump);
%     wl_seed = x_noPump(seedIdx);
%     seedPower(k, 1) = wl_seed;
%     seedPower(k, 2) = maxSeedPower;
% 
%     if abs(wl_seed - wl_pump) < 0.0015e-6
%         genPower(k,1) = NaN; genPower(k,2) = NaN;
%         continue;
%     end
% 
%     % Càlcul posició teòrica Idler (Generated Signal)
%     w_pump       = c / wl_pump;
%     w_seed       = c / wl_seed;
%     w_idler      = 2*w_pump - w_seed;
%     wl_idlerTheo = c / w_idler;
% 
%     int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window)) ...
%         & (abs(x - wl_pump) > pump_window) ...   
%         & (abs(x - wl_seed) > pump_window);      
% 
%     if any(int)
%         [maxY, relIdx] = max(y(int));
% 
%         lo_band = (x >= (wl_idlerTheo - 3*search_window)) & (x < (wl_idlerTheo - search_window)) ...
%                 & (abs(x - wl_pump) > pump_window) & (y > OSA_floor);
%         hi_band = (x >  (wl_idlerTheo + search_window)) & (x <= (wl_idlerTheo + 3*search_window)) ...
%                 & (abs(x - wl_pump) > pump_window) & (y > OSA_floor);
%         flank_y = y(lo_band | hi_band);
% 
%         if numel(flank_y) >= 5
%             local_noise_floor = median(flank_y);
%         else
%             valid_int = y(int) > -200;
%             if any(valid_int)
%                 local_noise_floor = median(y(int(valid_int)));
%             else
%                 local_noise_floor = -80; 
%             end
%         end
% 
%         if (maxY - local_noise_floor) >= snr_threshold
%             idx = find(int);
%             genPower(k, 1) = x(idx(relIdx));
%             genPower(k, 2) = maxY;
%         else
%             genPower(k, 1) = NaN; genPower(k, 2) = NaN;
%         end
%     else
%         genPower(k, 1) = NaN; genPower(k, 2) = NaN;
%     end
% end
% 
% % =========================================================================
% %% SECCIÓ 4: NETEJA I ORDENACIÓ DE DADES
% % =========================================================================
% valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);
% wl_idler_valid = genPower(valid, 1) * 1e9;
% P_idler_dBm    = genPower(valid, 2);
% P_seed_dBm     = seedPower(valid, 2);
% 
% [wl_idler_valid, sortIdx] = sort(wl_idler_valid);
% P_idler_dBm = P_idler_dBm(sortIdx);
% 
% % Neteja paral·lela de la llavor
% seed_wl_nm = seedPower(:,1) * 1e9;
% seed_pwr = seedPower(:,2);
% valid_seed_idx = abs(seed_wl_nm - 1550.12) > 0.5 & ~isnan(seed_wl_nm);
% clean_seed_wl_nm = seed_wl_nm(valid_seed_idx);
% clean_seed_pwr = seed_pwr(valid_seed_idx);
% 
% % =========================================================================
% %% SECCIÓ 5: FIGURA PÒSTER - SUBPLOTS SEED I GENERATED SIGNAL (Pastel)
% % =========================================================================
% fig_subplots = figure('Name', 'Power Profiling Comparison', 'Position', [100, 100, 800, 650]);
% 
% % --- SUBPLOT 1: Input Seed Power Sweep (Pastel Blue) ---
% subplot(2, 1, 1);
% plot(clean_seed_wl_nm, clean_seed_pwr, '-s', 'LineWidth', 1.5, ...
%      'Color', '#90CAF9', 'MarkerFaceColor', '#90CAF9', 'MarkerSize', 4);
% grid on;
% xlabel('Seed Wavelength (nm)', 'FontSize', 10);
% ylabel('Seed Power (dBm)', 'FontSize', 10);
% title('Input Seed Power Sweep (Fabry-Perot)', 'FontSize', 11);
% 
% % --- SUBPLOT 2: Generated Signal Power vs Wavelength (Pastel Red) ---
% subplot(2, 1, 2);
% plot(wl_idler_valid, P_idler_dBm, '-s', 'LineWidth', 1.5, ...
%      'Color', '#EF9A9A', 'MarkerFaceColor', '#EF9A9A', 'MarkerSize', 4);
% grid on;
% xlabel('Generated Photon Wavelength (nm)', 'FontSize', 10);
% ylabel('Generated Photon Power (dBm)', 'FontSize', 10);
% title('Generated Signal Power Performance ($wvg_2$ TE)', 'FontSize', 11);
% 
% % Desem exclusivament aquesta imatge unificada
% saveas(fig_subplots, fullfile(outFolder, 'Seed_vs_Generated_Subplots.png'));
% fprintf('\nProcés completat! S''ha generat el gràfic compost a: %s\n', fullfile(outFolder, 'Seed_vs_Generated_Subplots.png'));








































%%

% =========================================================================
% TFG - EVOLUCIÓ DE L'ESPECTRE FWM (Estil Publicació Científica / Pòster)
% =========================================================================
% ÍNDEX DEL CODI:
%   SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
%   SECCIÓ 2: CÀRREGA DE DADES (Load Data de wvg2_TM1)
%   SECCIÓ 3: SQUEEZE PIPELINE (Peak Tracking mínim necessari per a màscares)
%   SECCIÓ 4: GENERACIÓ DEL PLOT OVERLAY PROFESSIONAL (Blue-to-Magenta)
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
%% SECCIÓ 2: CÀRREGA DE DADES (Load Data)
% =========================================================================
dataFolder = '../Tests/stimFWM/wvg2_TE1';
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);
sweepData = cell(numFiles, 1);

fprintf('Carregant %d fitxers de la carpeta %s...\n', numFiles, dataFolder);
for k = 1:numFiles
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    sweepData{k} = readmatrix(fullFileName);
end

outFolder = dataFolder; 
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

% =========================================================================
%% SECCIÓ 3: SQUEEZE PIPELINE (Mínim Tracking per localitzar el Seed)
% =========================================================================
% Constants necessàries per fer la neteja de traces
wl_pump = 1.5496e-06;
pump_window = 0.5e-9;
seedPower = zeros(numFiles, 2);

fprintf('Executant rastreig ràpid de posicions de llavor (Seed Tracking)...\n');
for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);
    if mean(x) > 1000, x = x * 1e-9; end
    
    % Cerca de la llavor fora de la zona del pump per crear les màscares després
    not_pump  = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
    x_noPump  = x(not_pump);
    y_noPump  = y(not_pump);
    
    if ~isempty(y_noPump)
        [maxSeedPower, seedIdx] = max(y_noPump);
        seedPower(k, 1) = x_noPump(seedIdx);
        seedPower(k, 2) = maxSeedPower;
    else
        seedPower(k, 1) = NaN;
        seedPower(k, 2) = NaN;
    end
end

% =========================================================================
%% SECCIÓ 4: GENERACIÓ DEL PLOT OVERLAY PROFESSIONAL (Blue-to-Magenta)
% =========================================================================
fprintf('Generant el gràfic d''evolució espectral optimitzat...\n');
fig_overlay_clean = figure('Name', 'Spectra Overlay Publication', 'Position', [150, 150, 950, 520]);
hold on;

% Subconjunt indexat de traces a dibuixar (les de la teva òrbita/escombrat)
subset_files = [1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 20, 30, 35, 40, 50, 70, 72, 75, 80, 85, 90, 95, 100, 110, 120, 130, 135, 140, 146];
nLines = length(subset_files);

% Crear un gradient binari altament professional: Deep Blue a Vibrant Magenta
c1 = [0.12, 0.47, 0.71]; % Blau fosc inicial
c2 = [0.73, 0.15, 0.48]; % Magenta/Lila final
custom_map = [linspace(c1(1),c2(1),nLines)', linspace(c1(2),c2(2),nLines)', linspace(c1(3),c2(3),nLines)'];

for j = 1:nLines
    idx = subset_files(j);
    if idx <= numFiles
        x_nm = sweepData{idx}(:,1);
        y_dBm = sweepData{idx}(:,2);
        if mean(x_nm) < 1000, x_nm = x_nm * 1e9; end
        
        % Aplicació de la màscara per "tallar" el pic del Seed d'aquesta traça
        wl_seed_actual_nm = seedPower(idx, 1) * 1e9;
        if ~isnan(wl_seed_actual_nm)
            seed_mask = (x_nm >= (wl_seed_actual_nm - 1.5)) & (x_nm <= (wl_seed_actual_nm + 1.5));
            y_dBm(seed_mask) = NaN;
        end
        
        % Traçat net amb l'amplada de banda òptima per a pòsters
        plot(x_nm, y_dBm, 'Color', custom_map(j,:), 'LineWidth', 1.3);
    end
end

% Referència fina i elegant del Pump central (estil puntejat discret)
xline(wl_pump*1e9, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8, ...
      'Label', 'Pump', 'LabelOrientation', 'horizontal', 'FontSize', 10);

% Configuració d'etiquetes en LaTeX
xlabel('Wavelength (nm)', 'FontSize', 12); 
ylabel('Power (dBm)', 'FontSize', 12);
title('Spectral Evolution of the Generated Signal', 'FontSize', 13, 'FontWeight', 'bold');

% Estilització externa de la caixa de l'eix (Normes d'impacte científic)
box on;
set(gca, ...
    'TickDir', 'out', ...          % Ticks apuntant cap a fora
    'LineWidth', 1.2, ...          % Marc sòlid
    'FontSize', 11, ...            % Text altament llegible per al pòster
    'GridColor', [0.85 0.85 0.85], ... % Reixeta atenuada de fons
    'GridAlpha', 0.5);
grid on; 
ylim([-90 -20]);

% Implementació de la barra de color lateral per indicar evolució temporal
colormap(custom_map);
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Experimental Sweep Index ($\rightarrow$ Time Evolution)';
cb.Label.FontSize = 11;
% Marcatge lògic basat en el teu vector real de fitxers
set(cb, 'Ticks', [0, 0.5, 1], 'TickLabels', {num2str(subset_files(1)), num2str(median(subset_files)), num2str(subset_files(end))});

hold off;

% Desat en alta resolució PNG a la mateixa carpeta de dades
output_name = fullfile(outFolder, 'Spectra_Overlay_Publication.png');
saveas(fig_overlay_clean, output_name);
fprintf('\n[ÈXIT] S''ha completat l''execució en temps rècord.\nImatge desada a: %s\n', output_name);