clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%-------------
%  LOAD DATA
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

%%
% ------------
% SINGLE PLOT
% ------------
sanityCheck = sweepData{1};
figure(1)
plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');

sanityCheck = sweepData{11};
figure(2)
plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000')

sanityCheck = sweepData{31};
figure(3)
plot(sanityCheck(:,1), sanityCheck(:,2), 'Color', '#000000');