%% SECURITY 1
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
dataFolder = '../Tests/Chip13_wvg1_TEstimFWM_test_5';
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles=length(csvFiles);

sweepData = cell(numFiles, 1);
for k = 1:numFiles
    % Path for current file
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    % Load data
    sweepData{k} = readmatrix(fullFileName);
end

% ------------
%% SINGLE PLOT
% ------------
%for i = [1, 11, 25, 30, 50, 51, 60, 62, 65, 70, 80, 85, 90, 98, 99]
%for i = [1, 2, 3, 11, 16, 20, 50, 60, 70, 80, 90, 100, 105, 110, 111]
for i = [1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 20, 30, 35, 40, 50, 70, 72, 75, 80, 85, 90, 95, 100, 110, 120, 130, 135, 140, 146]
%for i = [1, 2, 3, 10, 20, 40, 50, 60, 70, 71]
    sanityCheck = sweepData{i};
    fig = figure(i);
    plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');
    xlabel('Wavelength'); ylabel('Power (dBm)');
    title(sprintf('Sweep Data for File %d', i));
    grid on;
    outFolder = '../Tests/Chip13_wvg1_TEstimFWM_test_5';
    % if ~exist(outFolder, 'dir'), mkdir(outFolder); end
    % fname = fullfile(outFolder, sprintf('plot_single_%02d.png', i));
    % saveas(fig, fname);
    i = i+1;
end

% ------------------------
%% PEAK TRACKING PIPELINE
% ------------------------
% Physical constants
c = 299792458;
wl_pump = 1.5496e-06; 
pump_window = 0.5e-9;       
search_window = 1.0e-9;     

% --- Noise threshold---
snr_threshold = 133;        
% --------------------------------

genPower = zeros(numFiles, 2);
seedPower = zeros(numFiles, 2);
pumpPower = zeros(numFiles, 1); % --- NEW: Vector to store Pump Power ---

for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);
    
    if mean(x) > 1000
        x = x * 1e-9;
    end
    
    % --- NEW: Extract pump power ---
    is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
    if any(is_pump)
        pumpPower(k) = max(y(is_pump));
    else
        pumpPower(k) = NaN;
    end
    % -------------------------------

    % Find seed position
    not_pump = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
    x_noPump = x(not_pump);
    y_noPump = y(not_pump);
    
    [maxSeedPower, seedIdx] = max(y_noPump);
    wl_seed = x_noPump(seedIdx);
    seedPower(k, 1) = wl_seed; 
    seedPower(k, 2) = maxSeedPower; 
    
    % Reject seed = pump
    if abs(wl_seed - wl_pump) < 0.0015e-6
        genPower(k,1) = NaN;
        genPower(k,2) = NaN;
        continue;           
    end
    
    % Calculate idler theoretical position
    w_pump = c / wl_pump; 
    w_seed = c / wl_seed;
    w_idler = 2*w_pump - w_seed;
    wl_idlerTheo = c / w_idler;
    
    int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window));
    
    if any(int)
        [maxY, relIdx] = max(y(int));

        lo_band = (x >= (wl_idlerTheo - 3*search_window)) & (x < (wl_idlerTheo - search_window));
        hi_band = (x >  (wl_idlerTheo + search_window))   & (x <= (wl_idlerTheo + 3*search_window));
        flank_y = y(lo_band | hi_band);
        if numel(flank_y) >= 5
            local_noise_floor = median(flank_y);
        % else
        %     local_noise_floor = median(y(int));  % fallback if flanks are outside OSA range
        end
        
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
disp(snr_threshold);
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
% 1. Clean up the data (Robust check: require valid Idler AND Pump AND Seed)
valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);

wl_idler_valid = genPower(valid, 1) * 1e9;
P_idler_dBm = genPower(valid, 2);
P_seed_dBm = seedPower(valid, 2);
P_pump_dBm = pumpPower(valid);

% Sort arrays by wavelength to prevent zig-zag plots
[wl_idler_valid, sortIdx] = sort(wl_idler_valid);
P_idler_dBm = P_idler_dBm(sortIdx);
P_seed_dBm  = P_seed_dBm(sortIdx);
P_pump_dBm  = P_pump_dBm(sortIdx);

% 2. Convert dBm to Watts (Linear scale)
P_idler_W = 10.^((P_idler_dBm - 30) / 10);
P_seed_W  = 10.^((P_seed_dBm - 30) / 10);
P_pump_W  = 10.^((P_pump_dBm - 30) / 10);

% 3. Calculate Normalized Efficiency (W^-2)
norm_eff_W2 = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));

% 4. Convert back to log scale for plotting
norm_eff_dB_W2 = 10 * log10(norm_eff_W2);

% 5. Plot
fig = figure('Name', 'Normalized FWM Efficiency');
plot(wl_idler_valid, norm_eff_dB_W2, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
grid on;
xlabel('Generated photon wavelength (nm)'); 
ylabel('Normalized Efficiency $10\log_{10}(\eta_{norm})$ [dB W$^{-2}$]');
title('Normalized FWM Efficiency ($\eta_{norm}$)');
ylim([0 100])

% 6. Save
outFolder = '../Tests/Chip13_wvg1_TEstimFWM_test_5';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig, fullfile(outFolder, 'Norm_ce.png'));

















%% SECURITY 2
% --------------------------------------------------------------
%% Power vs wvl (from OSA)
%% FWM Normalized Conversion Efficiency
% --------------------------------------------------------------

% figure(198); A=sweepData{98}; plot(A(:,1), A(:,2))

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
dataFolder = '../Tests/Chip13_wvg2_TEstimFWM_test_4';
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles=length(csvFiles);

sweepData = cell(numFiles, 1);
for k = 1:numFiles
    % Path for current file
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    % Load data
    sweepData{k} = readmatrix(fullFileName);
end

% ------------
%% SINGLE PLOT
% ------------
%for i = [1, 11, 25, 30, 50, 51, 60, 62, 65, 70, 80, 85, 90, 98, 99]
%for i = [1, 2, 3, 11, 16, 20, 50, 60, 70, 80, 90, 100, 105, 110, 111]
for i = [1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 20, 30, 35, 40, 50, 70, 72, 75, 80, 85, 90, 95, 100, 110, 120, 130, 135, 140, 146]
%for i = [1, 2, 3, 10, 20, 40, 50, 60, 70, 71]
    sanityCheck = sweepData{i};
    fig = figure(i);
    plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');
    xlabel('Wavelength'); ylabel('Power (dBm)');
    title(sprintf('Sweep Data for File %d', i));
    grid on;
    %outFolder = '../Tests/Chip13_wvg2_TEstimFWM_test_4';
    % if ~exist(outFolder, 'dir'), mkdir(outFolder); end
    % fname = fullfile(outFolder, sprintf('plot_single_%02d.png', i));
    % saveas(fig, fname);
    i = i+1;
end

% ------------------------
%% PEAK TRACKING PIPELINE
% ------------------------
% Physical constants
c = 299792458;
wl_pump = 1.5496e-06; 
pump_window = 0.5e-9;       
search_window = 1.0e-9;     

% --- Noise threshold---
snr_threshold = 20;        
% --------------------------------

genPower = zeros(numFiles, 2);
seedPower = zeros(numFiles, 2);
pumpPower = zeros(numFiles, 1); % --- NEW: Vector to store Pump Power ---

for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);
    
    if mean(x) > 1000
        x = x * 1e-9;
    end
    
    % --- NEW: Extract pump power ---
    is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
    if any(is_pump)
        pumpPower(k) = max(y(is_pump));
    else
        pumpPower(k) = NaN;
    end
    % -------------------------------

    % Find seed position
    not_pump = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
    x_noPump = x(not_pump);
    y_noPump = y(not_pump);
    
    [maxSeedPower, seedIdx] = max(y_noPump);
    wl_seed = x_noPump(seedIdx);
    seedPower(k, 1) = wl_seed; 
    seedPower(k, 2) = maxSeedPower; 
    
    % Reject seed = pump
    if abs(wl_seed - wl_pump) < 0.0015e-6
        genPower(k,1) = NaN;
        genPower(k,2) = NaN;
        continue;           
    end
    
    % Calculate idler theoretical position
    w_pump = c / wl_pump; 
    w_seed = c / wl_seed;
    w_idler = 2*w_pump - w_seed;
    wl_idlerTheo = c / w_idler;
    
    int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window));
    
    if any(int)
        [maxY, relIdx] = max(y(int));

        lo_band = (x >= (wl_idlerTheo - 3*search_window)) & (x < (wl_idlerTheo - search_window)) ...
                & (abs(x - wl_pump) > pump_window);   % <-- exclou cua del pump
        hi_band = (x >  (wl_idlerTheo + search_window))   & (x <= (wl_idlerTheo + 3*search_window)) ...
                & (abs(x - wl_pump) > pump_window);   % <-- exclou cua del pump
        flank_y = y(lo_band | hi_band);
        if numel(flank_y) >= 5
            local_noise_floor = median(flank_y);
        else
            local_noise_floor = median(y(int));  % fallback
        end
        
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
disp(snr_threshold);
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
% 1. Clean up the data (Robust check: require valid Idler AND Pump AND Seed)
valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1)) & ~isnan(pumpPower);

wl_idler_valid = genPower(valid, 1) * 1e9;
P_idler_dBm = genPower(valid, 2);
P_seed_dBm = seedPower(valid, 2);
P_pump_dBm = pumpPower(valid);

% Sort arrays by wavelength to prevent zig-zag plots
[wl_idler_valid, sortIdx] = sort(wl_idler_valid);
P_idler_dBm = P_idler_dBm(sortIdx);
P_seed_dBm  = P_seed_dBm(sortIdx);
P_pump_dBm  = P_pump_dBm(sortIdx);

% 2. Convert dBm to Watts (Linear scale)
P_idler_W = 10.^((P_idler_dBm - 30) / 10);
P_seed_W  = 10.^((P_seed_dBm - 30) / 10);
P_pump_W  = 10.^((P_pump_dBm - 30) / 10);

% 3. Calculate Normalized Efficiency (W^-2)
norm_eff_W2 = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));

% 4. Convert back to log scale for plotting
norm_eff_dB_W2 = 10 * log10(norm_eff_W2);

% 5. Plot
fig = figure('Name', 'Normalized FWM Efficiency');
plot(wl_idler_valid, norm_eff_dB_W2, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
grid on;
xlabel('Generated photon wavelength (nm)'); 
ylabel('Normalized Efficiency $10\log_{10}(\eta_{norm})$ [dB W$^{-2}$]');
title('Normalized FWM Efficiency ($\eta_{norm}$)');
ylim([0 100])

% 6. Save
outFolder = '../Tests/Chip13_wvg2_TEstimFWM_test_4';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig, fullfile(outFolder, 'Norm_ce.png'));