% ----------
%% SETTINGS
% ----------
clear; close all; format long g
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


%% 4 waveguides different length
power = [8.2, 4.8, 1.45, -2.53];
length = [8800, 26235.508836, 43671.01763, 61106.526509] / 10000;

fig1 = figure('Name', 'Prop loss', 'Position', [100, 100, 700, 500]);
plot(length, power, 'xk', 'LineWidth', 1.5, 'MarkerSize',8);
xlabel('Waveguide Length (cm)'); ylabel('Power (dBm)'); grid on;
title('Power vs. Waveguide Length');

% Linear fit
p = polyfit(length, power, 1);
yfit = polyval(p, length);
hold on; plot(length, yfit, '-', 'LineWidth', 1, 'Color', '#a9a9a9');

p = polyfit(length, power, 1);  % p(1)=m, p(2)=n
m = p(1);
n = p(2);
disp('Linear fit values for propagation and coupling loss calculation:')
fprintf('m = %.6g\nn = %.6g\n', m, n);
disp()


% Save results
outFolder = '../Plots/Chip1';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig1, fullfile(outFolder, 'propagation_loss.png'));
