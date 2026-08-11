% =========================================================================
%%                     ENTANGLEMENT INTERFERENCE
% =========================================================================
clear all; close all; format long g; clc;

set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');


% *******************************
%%       Data processing
% *******************************
% DATA
tau = 400 * 1e-12; % coincidence window 400ps
data_path = '../../Tests/timetagger/wvg2/entr9';
files = dir(fullfile(data_path, '*A*-B*.txt')); %entr4
%files = dir(fullfile(data_path, '*INTERFERENCE_A-*_B-*.txt')); % entr2
%SWITCH ALSO ALICE AND BOB 12--13

% ARRAYS to hold data
ang_H = []; raw_HH_H = []; raw_HV_H = []; raw_VV_H = []; raw_VH_H = []; % HH = CC1; HV = CC3
ang_D = []; raw_HH_D = []; raw_HV_D = []; raw_VV_D = []; raw_VH_D = [];
ch1_H = []; ch2_H = []; ch3_H = []; ch4_H = []; % single counts
ch1_D = []; ch2_D = []; ch3_D = []; ch4_D = [];

% Load and process each file
for k = 1:length(files)
    file_path = fullfile(data_path, files(k).name);
    data = readmatrix(file_path, 'NumHeaderLines', 12);  % Read numeric data starting from line 13

    RM_Alice = data(1, 13); RM_Bob = data(1, 12)*2;

    if RM_Alice == 0.0 % H base
        % Single counts
        ch1_H = [ch1_H, data(1, 3)]; %HA*
        ch2_H = [ch2_H, data(1, 4)]; %VA
        ch3_H = [ch3_H, data(1, 5)]; %HB*
        ch4_H = [ch4_H, data(1, 6)]; %VB*

        % CC
        raw_HH_H = [raw_HH_H, data(1, 8)];
        raw_HV_H = [raw_HV_H, data(1, 10)];
        raw_VV_H = [raw_VV_H, data(1, 9)];
        raw_VH_H = [raw_VH_H, data(1, 11)];

        % Bob angle
        ang_H = [ang_H, data(1, 12)*2];

        % Accidentals calculation
        acc_HH = ch1_H .* ch3_H * tau;
        acc_HV = ch1_H .* ch4_H * tau;
        acc_VV = ch2_H .* ch4_H * tau;
        acc_VH = ch2_H .* ch3_H * tau;

        % Real coincidences
        HH_H = raw_HH_H - acc_HH;
        HV_H = raw_HV_H - acc_HV;
        VV_H = raw_VV_H - acc_VV;
        VH_H = raw_VH_H - acc_VH;
    end
    if RM_Alice== 22.5
        % Single counts
        ch1_D = [ch1_D, data(1, 3)]; %HA*
        ch2_D = [ch2_D, data(1, 4)]; %VA
        ch3_D = [ch3_D, data(1, 5)]; %HB*
        ch4_D = [ch4_D, data(1, 6)]; %VB*

        % CC
        raw_HH_D = [raw_HH_D, data(1, 8)];
        raw_HV_D = [raw_HV_D, data(1, 10)];
        raw_VV_D = [raw_VV_D, data(1, 9)];
        raw_VH_D = [raw_VH_D, data(1, 11)];

        % Bob angle
        ang_D = [ang_D, data(1, 12)*2];

        % Accidentals calculation
        acc_HH = ch1_D .* ch3_D * tau;
        acc_HV = ch1_D .* ch4_D * tau;
        acc_VV = ch2_D .* ch4_D * tau;
        acc_VH = ch2_D .* ch3_D * tau;

        % Real coincidences
        HH_D = raw_HH_D - acc_HH;
        HV_D = raw_HV_D - acc_HV;
        VV_D = raw_VV_D - acc_VV;
        VH_D = raw_VH_D - acc_VH;
    end
end

% Sort by angle in ascending order
[ang_H_sort, idxH] = sort(ang_H, 'ascend');
HH_H_sort = HH_H(idxH);
HV_H_sort = HV_H(idxH);

[ang_D_sort, idxD] = sort(ang_D, 'ascend');
HH_D_sort = HH_D(idxD);
HV_D_sort = HV_D(idxD);

% Calculate raw (max-min) visibility, kept for comparison with fitted value
vis_H = (max(HH_H_sort) - min(HH_H_sort)) / (max(HH_H_sort) + min(HH_H_sort));
vis_D = (max(HH_D_sort) - min(HH_D_sort)) / (max(HH_D_sort) + min(HH_D_sort));

fprintf('Alice at H --> Raw Visibility: %.2f%%\n', vis_H*100);
fprintf('Alice at D --> Raw Visibility: %.2f%%\n', vis_D*100);



% *******************************
%%       Cosine fit
% *******************************
% Model: y = a*cos(b*x + c) + d
% Physical prior: rotating the analyzer produces a fringe with period
% 180 deg in the Bob angle (Malus's law -> cos^2 -> period pi in angle*2).
% This fixes b0 sensibly instead of guessing blindly, which is what was
% causing lsqcurvefit to lock onto bad local minima before.

cosineModel = @(p, x) p(1) * cos(p(2) * x + p(3)) + p(4);
opts = optimoptions('lsqcurvefit', 'Display', 'off');

fitCosine = @(x, y) local_fit_cosine(x, y, cosineModel, opts);

[pFit_HH_H, xq_H, yq_HH_H] = fitCosine(ang_H_sort, HH_H_sort);
[pFit_HV_H, ~,    yq_HV_H] = fitCosine(ang_H_sort, HV_H_sort);
[pFit_HH_D, xq_D, yq_HH_D] = fitCosine(ang_D_sort, HH_D_sort);
[pFit_HV_D, ~,    yq_HV_D] = fitCosine(ang_D_sort, HV_D_sort);

% Fitted visibility from each HH fit: V = a / d
visFit_H = abs(pFit_HH_H(1)) / pFit_HH_H(4);
visFit_D = abs(pFit_HH_D(1)) / pFit_HH_D(4);

fprintf('Alice at H --> Fitted Visibility (HH): %.2f%%\n', visFit_H*100);
fprintf('Alice at D --> Fitted Visibility (HH): %.2f%%\n', visFit_D*100);



% ********************************************
%%       Plot 1: divided by Alice base
% ********************************************
figure('Name','Interference Fringes','Color','w','Position',[100 100 1200 500]);

col_HH = [0.35, 0.60, 0.90];   % pastel blue
col_HV = [0.92, 0.45, 0.45];   % pastel red

% -------- Plot 1: Alice at H --------
ax1 = subplot(1,2,1);
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

plot(ax1, xq_H, yq_HH_H, '-', 'Color', col_HH, 'LineWidth', 1.5, ...
    'DisplayName', 'HH fit');
plot(ax1, xq_H, yq_HV_H, '-', 'Color', col_HV, 'LineWidth', 1.5, ...
    'DisplayName', 'HV fit');

plot(ax1, ang_H_sort, HH_H_sort, 'o', ...
    'MarkerFaceColor', col_HH, 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'DisplayName', 'HH data');
plot(ax1, ang_H_sort, HV_H_sort, 's', ...
    'MarkerFaceColor', col_HV, 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'DisplayName', 'HV data');

xlabel(ax1, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax1, 'True Coincidences (cps)', 'FontSize', 12);
title(ax1, sprintf('Alice at $H$ | $V_{fit} = %.1f\\%%$', visFit_H*100), 'FontSize', 13);
legend(ax1, 'Location', 'best', 'FontSize', 10);
set(ax1, 'FontSize', 11);

% -------- Plot 2: Alice at D --------
ax2 = subplot(1,2,2);
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

plot(ax2, xq_D, yq_HH_D, '-', 'Color', col_HH, 'LineWidth', 1.5, ...
    'DisplayName', 'HH fit');
plot(ax2, xq_D, yq_HV_D, '-', 'Color', col_HV, 'LineWidth', 1.5, ...
    'DisplayName', 'HV fit');

plot(ax2, ang_D_sort, HH_D_sort, 'o', ...
    'MarkerFaceColor', col_HH, 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'DisplayName', 'HH data');
plot(ax2, ang_D_sort, HV_D_sort, 's', ...
    'MarkerFaceColor', col_HV, 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, 'DisplayName', 'HV data');

xlabel(ax2, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax2, 'True Coincidences (cps)', 'FontSize', 12);
title(ax2, sprintf('Alice at $D$ | $V_{fit} = %.1f\\%%$', visFit_D*100), 'FontSize', 13);
legend(ax2, 'Location', 'best', 'FontSize', 10);
set(ax2, 'FontSize', 11);

% Print all fitted parameters for reference/debugging
fprintf('\n--- Fit parameters [a, b, c, d] ---\n');
fprintf('HH_H: a=%.3f, b=%.5f rad/deg, c=%.3f rad, d=%.3f, T=%.2f deg\n', pFit_HH_H, 2*pi/pFit_HH_H(2));
fprintf('HV_H: a=%.3f, b=%.5f rad/deg, c=%.3f rad, d=%.3f, T=%.2f deg\n', pFit_HV_H, 2*pi/pFit_HV_H(2));
fprintf('HH_D: a=%.3f, b=%.5f rad/deg, c=%.3f rad, d=%.3f, T=%.2f deg\n', pFit_HH_D, 2*pi/pFit_HH_D(2));
fprintf('HV_D: a=%.3f, b=%.5f rad/deg, c=%.3f rad, d=%.3f, T=%.2f deg\n', pFit_HV_D, 2*pi/pFit_HV_D(2));




% ********************************************
%%       Plot 2: All together
% ********************************************
figure('Name','Interference Fringes','Color','w','Position',[100 100 1200 500]);

% Pastel colors (define separate shades for H vs D)
col_HH_H = [0.65, 0.85, 1.00];  % pastel blue for Alice at H
col_HH_D = [0.35, 0.60, 0.90];  % pastel blue for Alice at D

col_HV_H = [1.00, 0.65, 0.65];  % pastel red for Alice at H
col_HV_D = [0.92, 0.45, 0.45];  % pastel red for Alice at D

ax = axes; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

% ---------------- Fits (NOT in legend) ----------------
% HH fits
pHH_Hfit = plot(ax, xq_H, yq_HH_H, '-', 'Color', col_HH_H, 'LineWidth', 1.5, ...
    'HandleVisibility','off'); % hides from legend
pHH_Dfit = plot(ax, xq_D, yq_HH_D, '-', 'Color', col_HH_D, 'LineWidth', 1.5, ...
    'HandleVisibility','off');

% HV fits
pHV_Hfit = plot(ax, xq_H, yq_HV_H, '-', 'Color', col_HV_H, 'LineWidth', 1.5, ...
    'HandleVisibility','off');
pHV_Dfit = plot(ax, xq_D, yq_HV_D, '-', 'Color', col_HV_D, 'LineWidth', 1.5, ...
    'HandleVisibility','off');

% ---------------- Experimental data (IN legend) ----------------
plot(ax, ang_H_sort, HH_H_sort, 'o', ...
    'MarkerFaceColor', col_HH_H, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','HH (Alice at H)');

plot(ax, ang_H_sort, HV_H_sort, 's', ...
    'MarkerFaceColor', col_HV_H, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','HV (Alice at H)');

plot(ax, ang_D_sort, HH_D_sort, 'o', ...
    'MarkerFaceColor', col_HH_D, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','HH (Alice at D)');

plot(ax, ang_D_sort, HV_D_sort, 's', ...
    'MarkerFaceColor', col_HV_D, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','HV (Alice at D)');

% ---------------- Labels/title/legend ----------------
xlabel(ax, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax, 'True Coincidences (cps)', 'FontSize', 12);

title(ax, sprintf('Interference Fringes. $V_{fit}^H = %.1f\\%%$, $V_{fit}^D = %.1f\\%%$', ...
    visFit_H*100, visFit_D*100), 'FontSize', 13);

legend(ax, 'Location','best', 'FontSize', 10);
set(ax, 'FontSize', 11);




% ********************************************
%%       Plot 3: divided by HH / HV
% ********************************************
figure('Name','Interference Fringes','Color','w','Position',[100 100 1200 650]);

% Pastel colors (separate shades for H vs D)
col_HH_H = [0.65, 0.85, 1.00];  % pastel blue for Alice at H
col_HH_D = [0.35, 0.60, 0.90];  % pastel blue for Alice at D

col_HV_H = [1.00, 0.65, 0.65];  % pastel red for Alice at H
col_HV_D = [0.92, 0.45, 0.45];  % pastel red for Alice at D

% ---------- Top subplot: HH ----------
ax1 = subplot(2,1,1);
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

% Fits (NOT in legend)
plot(ax1, xq_H, yq_HH_H, '-', 'Color', col_HH_H, 'LineWidth', 1.5, ...
    'HandleVisibility','off');
plot(ax1, xq_D, yq_HH_D, '-', 'Color', col_HH_D, 'LineWidth', 1.5, ...
    'HandleVisibility','off');

% Experimental data (IN legend)
plot(ax1, ang_H_sort, HH_H_sort, 'o', ...
    'MarkerFaceColor', col_HH_H, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','Alice at H');

plot(ax1, ang_D_sort, HH_D_sort, 'o', ...
    'MarkerFaceColor', col_HH_D, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','Alice at D');

xlabel(ax1, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax1, 'True Coincidences (cps)', 'FontSize', 12);
title(ax1, sprintf('HH Coincidence Counts', ...
    visFit_H*100, visFit_D*100), 'FontSize', 13);

legend(ax1, 'Location','best', 'FontSize', 10);
set(ax1, 'FontSize', 11);

% ---------- Bottom subplot: HV ----------
ax2 = subplot(2,1,2);
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

% Fits (NOT in legend)
plot(ax2, xq_H, yq_HV_H, '-', 'Color', col_HV_H, 'LineWidth', 1.5, ...
    'HandleVisibility','off');
plot(ax2, xq_D, yq_HV_D, '-', 'Color', col_HV_D, 'LineWidth', 1.5, ...
    'HandleVisibility','off');

% Experimental data (IN legend)
plot(ax2, ang_H_sort, HV_H_sort, 's', ...
    'MarkerFaceColor', col_HV_H, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','Alice at H');

plot(ax2, ang_D_sort, HV_D_sort, 's', ...
    'MarkerFaceColor', col_HV_D, 'MarkerEdgeColor','k', ...
    'MarkerSize',7, 'DisplayName','Alice at D');

xlabel(ax2, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax2, 'True Coincidences (cps)', 'FontSize', 12);
title(ax2, sprintf('HV Coincidence Counts', ...
    visFit_H*100, visFit_D*100), 'FontSize', 13);

legend(ax2, 'Location','best', 'FontSize', 10);
set(ax2, 'FontSize', 11);



% =========================================================================
%% LOCAL FUNCTION
% =========================================================================
function [pFit, xq, yq] = local_fit_cosine(x, y, cosineModel, opts)
    % Build a physically-motivated initial guess and fit y = a*cos(b*x+c)+d

    x = x(:); y = y(:);

    d0 = mean(y);
    a0 = (max(y) - min(y)) / 2;
    b0 = 2*pi/180;                  % expected period = 180 deg (Malus's law)

    [~, idx_max] = max(y);
    c0 = -b0 * x(idx_max);

    p0 = [a0, b0, c0, d0];

    % Keep b close to the physically expected frequency to avoid
    % lsqcurvefit locking onto a spurious fast/slow oscillation.
    lb = [-Inf, 0.5*b0, -Inf, -Inf];
    ub = [ Inf, 1.5*b0,  Inf,  Inf];

    pFit = lsqcurvefit(cosineModel, p0, x, y, lb, ub, opts);

    xq = linspace(min(x), max(x), 300);
    yq = cosineModel(pFit, xq);
end