% --------------------------------------------------------------
%% Power vs wvl (from OSA)
%% FWM Conversion Efficiency
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
dataFolder = '../Tests/FWM_sweep_test_06';
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
for i = [1, 11, 30, 50, 51, 60, 80, 90, 100]
    sanityCheck = sweepData{i};
    fig = figure(i);
    plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');
    xlabel('Wavelength'); ylabel('Power (dBm)');
    title(sprintf('Sweep Data for File %d', i));
    grid on;
    outFolder = '../Tests/FWM_sweep_test_06';
    if ~exist(outFolder, 'dir'), mkdir(outFolder); end
    fname = fullfile(outFolder, sprintf('plot_single_%02d.png', i));
    saveas(fig, fname);
    i = i+1;
end


% ------------------------
%% PEAK TRACKING PIPELINE
% ------------------------

% Physical constants
c = 299792458;
wl_pump = 1.5496e-06; % exp 1550.12nm, but checked graphs
pump_window = 0.1e-9; %1.0e-9
search_window = 0.5e-9;

genPower = zeros(numFiles, 2);
seedPower = zeros(numFiles, 2);

for k = 1:numFiles
    x = sweepData{k}(:,1);
    y = sweepData{k}(:,2);

    % Find seed position. We ignore pump region.
    not_pump = (x < (wl_pump - pump_window)) | (x > (wl_pump + pump_window));
    x_noPump = x(not_pump);
    y_noPump = y(not_pump);

    % Max oustide region = seed
    [maxSeedPower, seedIdx] = max(y_noPump);
    wl_seed = x_noPump(seedIdx);
    seedPower(k, 1) = wl_seed; % Corresponding wavelength for seed
    seedPower(k, 2) = maxSeedPower; % Max seed power

    % Reject seed = pump
    if abs(wl_seed - wl_pump) < 8.8e-10 %0.001e-6%0.00112e-6%8.8e-10     % 8.8e-10 exact
        genPower(k,1) = NaN;
        genPower(k,2) = NaN;
        continue;           % maybe include wvl?
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
        genPower(k, 1) = x(idx(relIdx)); % Corresponding wavelength for generated power
        genPower(k, 2) = maxY; % Max generated power
    else
        genPower(k, 1) = NaN;
        genPower(k, 2) = NaN;
    end
end
%genPower(11,:) = []; seedPower(11,:) = [];

% -------------------
%% FWM VISUALIZATION
% -------------------

% 1. Clean up the data (Filter out NaNs so MATLAB can plot properly)
valid = ~isnan(genPower(:,1)) & ~isnan(seedPower(:,1));

wl_idler_valid = genPower(valid, 1) * 1e9;
P_idler_dBm_valid = genPower(valid, 2);
P_seed_dBm_valid = seedPower(valid, 2);

% 2. Calculate Efficiency
conv_eff_dB = P_idler_dBm_valid - P_seed_dBm_valid;

% (Optional linear calculation if you need it later)
genPower_mW = 10.^(P_idler_dBm_valid / 10);
seedPower_mW = 10.^(P_seed_dBm_valid / 10);
conv_eff_lin = genPower_mW ./ seedPower_mW;

% 3. Plot
fig = figure('Name', 'FWM Conversion Efficiency dB');
plot(wl_idler_valid, conv_eff_dB, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
grid on;
xlabel('Generated photon wavelength (nm)'); 
ylabel('FWM Conversion Efficiency (dB)')
title('FWM Conversion Efficiency')

% 4. Save (Make sure it saves to Test_06 folder, not Test_02)
outFolder = '../Tests/FWM_sweep_test_06';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig, fullfile(outFolder, 'FWM_conv_eff_dB.png'));