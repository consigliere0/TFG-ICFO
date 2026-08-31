%
% =========================================================================
% TFG - AMITJÀ DE TESI: SUBPLOTS DE POTÈNCIA (Waveguide 2 TE)
% =========================================================================
% ÍNDEX DEL CODI:
%   SECCIÓ 1: CONFIGURACIÓ INICIAL (Settings)
%   SECCIÓ 2: CÀRREGA DE DADES (Load Data)
%   SECCIÓ 3: EXTRACCIÓ DE PICS (Peak Tracking Pipeline)
%   SECCIÓ 4: NETEJA I ORDENACIÓ DE DADES
%   SECCIÓ 5: FIGURA PÒSTER - SUBPLOTS SEED I GENERATED SIGNAL (Pastel)
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
dataFolder = '../../Tests/stimFWM/wvg2_TE1';
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);
sweepData = cell(numFiles, 1);

for k = 1:numFiles
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    sweepData{k} = readmatrix(fullFileName);
end

outFolder = dataFolder; 
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

% =========================================================================
%% SECCIÓ 3: EXTRACCIÓ DE PICS (Peak Tracking Pipeline)
% =========================================================================
c = 299792458;
wl_pump = 1.5496e-06;
pump_window = 0.5e-9;
search_window = 1.0e-9;
snr_threshold = 3.2;        
OSA_floor     = -80;        

genPower  = zeros(numFiles, 2);
seedPower = zeros(numFiles, 2);
pumpPower = zeros(numFiles, 1);

for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);

    if mean(x) > 1000, x = x * 1e-9; end

    % Extreure Pump
    is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
    if any(is_pump)
        pumpPower(k) = max(y(is_pump));
    else
        pumpPower(k) = NaN;
    end

    % Extreure Seed
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

    % Càlcul posició teòrica Idler (Generated Signal)
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

% =========================================================================
%% SECCIÓ 4: NETEJA I ORDENACIÓ DE DADES
% =========================================================================
valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);
wl_idler_valid = genPower(valid, 1) * 1e9;
P_idler_dBm    = genPower(valid, 2);
P_seed_dBm     = seedPower(valid, 2);

[wl_idler_valid, sortIdx] = sort(wl_idler_valid);
P_idler_dBm = P_idler_dBm(sortIdx);

% Neteja paral·lela de la llavor
seed_wl_nm = seedPower(:,1) * 1e9;
seed_pwr = seedPower(:,2);
valid_seed_idx = abs(seed_wl_nm - 1550.12) > 0.5 & ~isnan(seed_wl_nm);
clean_seed_wl_nm = seed_wl_nm(valid_seed_idx);
clean_seed_pwr = seed_pwr(valid_seed_idx);


%% Seed vs wavelength & generated signal
figure(1)
plot(clean_seed_wl_nm, clean_seed_pwr, '-s', 'LineWidth', 1.5, ...
     'Color', '#90CAF9', 'MarkerFaceColor', '#90CAF9', 'MarkerSize', 4);
grid on;
xlabel('Seed Wavelength (nm)', 'FontSize', 10);
ylabel('Seed Off-Chip Power (dBm)', 'FontSize', 10);


figure(2)
plot(wl_idler_valid, P_idler_dBm, '-s', 'LineWidth', 1.5, ...
     'Color', '#EF9A9A', 'MarkerFaceColor', '#EF9A9A', 'MarkerSize', 4);
grid on;
xlabel('Generated Photon Wavelength (nm)', 'FontSize', 10);
ylabel('Generated Photon Off-Chip Power (dBm)', 'FontSize', 10); 
