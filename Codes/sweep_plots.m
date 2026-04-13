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
dataFolder = '../Tests/FWM_Sweep_Test_01';
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
for i = [1, 10, 11, 12, 31]
    sanityCheck = sweepData{i};
    figure(i)
    plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');
    xlabel('Wavelength'); ylabel('Power (dBm)')
    i = i+1;
end


% ------------------------
%% PEAK TRACKING PIPELINE
% ------------------------

% Physical constants
c = 299792458;
wl_pump = 1.54965e-06; % exp 1550.12nm, but checked graphs
pump_window = 1.0e-9;
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
    if abs(wl_seed - wl_pump) < 0.00112e-6%8.8e-10     % 8.8e-10 exact
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
genPower(11,:) = []; seedPower(11,:) = [];

% -------------------
%% FWM VISUALIZATION
% -------------------
wl_idler = genPower(:,1) * 1e9;

conv_eff_dB = genPower(:,2) - seedPower(:,2);

genPower_mW = 10.^(genPower(:,2) / 10);
seedPower_mW = 10.^(seedPower(:,2) / 10);
conv_eff_lin = genPower_mW ./ seedPower_mW;

figure('Name', 'FWM Conversion Efficiency dB')
plot(wl_idler, conv_eff_dB, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
axis([1530 1561 -40 0]); grid on
xlabel('Generated photon wavelength (nm)')
ylabel('FWM Conversion Efficiency (dB)')
title('FWM Conversion Efficiency')
