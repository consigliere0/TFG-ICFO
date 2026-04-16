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
dataFolder = '../Tests/FWM_sweep_test_02';      % Test 01 irrelevant
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

% ------------------------
%% PEAK TRACKING PIPELINE
% ------------------------
c = 299792458;
wl_pump = 1.54965e-06; 
pump_window = 0.1e-9; 
search_window = 0.5e-9;

genPower = zeros(numFiles, 2);
seedPower = zeros(numFiles, 2);
pumpPower = zeros(numFiles, 1); % New pump array

for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);
    
    % Extract pump power
    is_pump = (x >= (wl_pump - pump_window)) & (x <= (wl_pump + pump_window));
    if any(is_pump)
        pumpPower(k) = max(y(is_pump));
    else
        pumpPower(k) = NaN;
    end
    
    % Find seed position. We ignore pump region.
    not_pump = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
    x_noPump = x(not_pump);
    y_noPump = y(not_pump);
    
    % Max oustide region = seed
    [maxSeedPower, seedIdx] = max(y_noPump);
    wl_seed = x_noPump(seedIdx);
    seedPower(k, 1) = wl_seed; 
    seedPower(k, 2) = maxSeedPower; 
    
    % Reject seed = pump
    if abs(wl_seed - wl_pump) < 0.00112e-6
        genPower(k,1) = NaN;
        genPower(k,2) = NaN;
        continue;           
    end
    
    % Calculate idler theoretical position
    w_pump = c / wl_pump; w_seed = c / wl_seed;
    w_idler = 2*w_pump - w_seed;
    wl_idlerTheo = c / w_idler;
    
    % Dynamic window arround theo value
    int = (x >= (wl_idlerTheo - search_window)) & (x <= (wl_idlerTheo + search_window));
    if any(int)
        [maxY, relIdx] = max(y(int));
        idx = find(int);
        genPower(k, 1) = x(idx(relIdx)); 
        genPower(k, 2) = maxY; 
    else
        genPower(k, 1) = NaN;
        genPower(k, 2) = NaN;
    end
end

% Eliminar l'element que vas determinar als teus tests anteriors
% genPower(11,:) = []; 
% seedPower(11,:) = [];
% pumpPower(11) = [];

% ----------------------
%% FILTERING & ORDERING
% ----------------------
% Filter out NaNs
valid = ~isnan(genPower(:,1)) & ~isnan(pumpPower);

wl_idler_nm = genPower(valid, 1) * 1e9;
P_idler_dBm = genPower(valid, 2);
P_seed_dBm  = seedPower(valid, 2);
P_pump_dBm  = pumpPower(valid);

% Order x-axis from - to + (so that it doesn't zig zag)
[wl_idler_nm, sortIdx] = sort(wl_idler_nm);
P_idler_dBm = P_idler_dBm(sortIdx);
P_seed_dBm  = P_seed_dBm(sortIdx);
P_pump_dBm  = P_pump_dBm(sortIdx);

% ----------------------------------------
%% 3-dB BANDWIDTH & CONVERSION EFFICIENCY
% ----------------------------------------
conv_eff_dB = P_idler_dBm - P_seed_dBm;

fig1 = figure('Name', '3-dB Bandwidth', 'Position', [100, 100, 700, 500]);
plot(wl_idler_nm, conv_eff_dB, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
hold on;

% Trobar el pic d'eficiència
max_eff = max(conv_eff_dB);
threshold_3dB = max_eff - 3;

% Dibuixar referències
yline(max_eff, '--r', 'Max Efficiency', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
yline(threshold_3dB, '--b', '-3 dB Threshold', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');

% Cercar amplada de banda
idx_3dB = conv_eff_dB >= threshold_3dB;
if any(idx_3dB)
    wl_valid = wl_idler_nm(idx_3dB);
    bw_nm = max(wl_valid) - min(wl_valid);
    
    txt = sprintf('3-dB BW $\\approx$ %.2f nm', bw_nm);
    text(mean(wl_valid), threshold_3dB - 1, txt, 'BackgroundColor', 'w', 'EdgeColor', 'b', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
end

grid on;
xlabel('Generated photon wavelength (nm)'); 
ylabel('Conversion Efficiency (dB)');
title('FWM Conversion Efficiency \& 3-dB Bandwidth');
hold off;

% ------------------------------
%% NORMALIZED EFFICIENCY (W^-2)
% ------------------------------
% dBm --> Linear (Watts!! not mW)
P_idler_W = 10.^((P_idler_dBm - 30) / 10);
P_seed_W  = 10.^((P_seed_dBm - 30) / 10);
P_pump_W  = 10.^((P_pump_dBm - 30) / 10);

% Normalize
norm_eff_W2 = P_idler_W ./ (P_seed_W .* (P_pump_W.^2));

% Final result in log scale
norm_eff_dB_W2 = 10 * log10(norm_eff_W2);

fig2 = figure('Name', 'Normalized Efficiency', 'Position', [850, 100, 700, 500]);
plot(wl_idler_nm, norm_eff_dB_W2, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
grid on;
xlabel('Generated photon wavelength (nm)'); 
ylabel('Normalized Efficiency $10\log_{10}(\eta_{norm})$ [dB W$^{-2}$]');
title('Normalized FWM Efficiency ($\eta_{norm}$)');

% Save results
outFolder = dataFolder;
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig1, fullfile(outFolder, 'Bandwidth_3dB.png'));
saveas(fig2, fullfile(outFolder, 'Norm_ce.png'));