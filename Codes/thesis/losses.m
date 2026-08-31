%% Propagation and coupling loss
clear; close all; format long g

set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

%% Input power at chip facet
Pin = 14.6;   % dBm

%% Chip 1
%L1 = [8800, 26235.508836, 43671.01763] / 10000;   % cm
%P1 = [9.28, 5.10, 1.45];                          % dBm

P1 = [8.2, 4.8, 1.45, -2.53];
L1 = [8800, 26235.508836, 43671.01763, 61106.526509] / 10000;

%% Chip 2
% Replace with the measured values
L2 =  [8800, 26235.508836, 43671.01763, 61106.526509] / 10000;                     % cm
P2 = [9.28, 5.10, 1.45, -2.345];                             % dBm

%% Convert measured output powers to insertion loss
% Insertion loss = Pin - Pout
IL1 = Pin - P1;
IL2 = Pin - P2;

%% Linear fits
fit1 = polyfit(L1, IL1, 1);
fit2 = polyfit(L2, IL2, 1);

m1 = fit1(1);
b1 = fit1(2);

m2 = fit2(1);
b2 = fit2(2);

%% Extract losses
alpha_p1 = m1;
alpha_c1 = b1/2;

alpha_p2 = m2;
alpha_c2 = b2/2;

%% Average losses across the two characterized chips
alpha_p_mean = mean([alpha_p1, alpha_p2]);
alpha_c_mean = mean([alpha_c1, alpha_c2]);

% Chip-to-chip spread
alpha_p_range = max([alpha_p1, alpha_p2]) - min([alpha_p1, alpha_p2]);
alpha_c_range = max([alpha_c1, alpha_c2]) - min([alpha_c1, alpha_c2]);

%% Display results
fprintf('\nPropagation and coupling losses\n');
fprintf('--------------------------------\n');

fprintf('Chip 1:\n');
fprintf('  Propagation loss = %.4f dB/cm\n', alpha_p1);
fprintf('  Coupling loss    = %.4f dB/facet\n\n', alpha_c1);

fprintf('Chip 2:\n');
fprintf('  Propagation loss = %.4f dB/cm\n', alpha_p2);
fprintf('  Coupling loss    = %.4f dB/facet\n\n', alpha_c2);

fprintf('Representative value (mean of chips 1 and 2):\n');
fprintf('  Propagation loss = %.4f dB/cm\n', alpha_p_mean);
fprintf('  Coupling loss    = %.4f dB/facet\n\n', alpha_c_mean);

fprintf('Chip-to-chip range:\n');
fprintf('  Propagation loss = %.4f dB/cm\n', alpha_p_range);
fprintf('  Coupling loss    = %.4f dB/facet\n', alpha_c_range);

%% Figure
fig = figure('Name', 'Propagation and Coupling Loss', ...
             'Position', [100, 100, 750, 520]);

hold on;
grid on;
box on;

% Neutral colours
c1 = [0.25 0.25 0.25];
c2 = [0.45 0.55 0.60];

% Lighter colours for fitted lines
fitLightFactor = 0.45;
c1_fit = c1 + fitLightFactor * (1 - c1);
c2_fit = c2 + fitLightFactor * (1 - c2);

% Experimental points
h1 = plot(L1, IL1, 'o', ...
    'Color', c1, ...
    'MarkerFaceColor', c1, ...
    'MarkerSize', 7, ...
    'LineWidth', 1.2);

h2 = plot(L2, IL2, 's', ...
    'Color', c2, ...
    'MarkerFaceColor', c2, ...
    'MarkerSize', 7, ...
    'LineWidth', 1.2);

% Fit segments spanning only the measured data range
Lfit1 = linspace(min(L1), max(L1), 100);
Lfit2 = linspace(min(L2), max(L2), 100);

hfit1 = plot(Lfit1, polyval(fit1, Lfit1), '-', ...
    'Color', c1_fit, ...
    'LineWidth', 1.5);

hfit2 = plot(Lfit2, polyval(fit2, Lfit2), '-', ...
    'Color', c2_fit, ...
    'LineWidth', 1.5);

xlabel('Waveguide length, $L$ (cm)');
ylabel('Insertion Loss (dB)');

legend([h1, h2, hfit1, hfit2], ...
       {'Chip 1', 'Chip 2', ...
        sprintf('Chip 1 fit ($\\alpha_p=%.2f$ dB/cm)', alpha_p1), ...
        sprintf('Chip 2 fit ($\\alpha_p=%.2f$ dB/cm)', alpha_p2)}, ...
        'Location', 'southeast');

set(gca, 'FontSize', 12, 'LineWidth', 1);

% %% Save figure
% outFolder = '../../Plots/ChipComparison';
% 
% if ~exist(outFolder, 'dir')
%     mkdir(outFolder);
% end
% 
% exportgraphics(fig, ...
%     fullfile(outFolder, 'propagation_coupling_loss.png'), ...
%     'Resolution', 300);