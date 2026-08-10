% =========================================================================
%%                     ENTANGLEMENT INTERFERENCE 
% =========================================================================
clear all; close all; format long g;

set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

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
%HH_H_sort = HH_H(idxH);
HH_H_sort = HH_H(idxH);
HV_H_sort = HV_H(idxH);

[ang_D_sort, idxD] = sort(ang_D, 'ascend');
HH_D_sort = HH_D(idxD);
HV_D_sort = HV_D(idxD);


% Calculate visibility
vis_H = (max(HH_H_sort) - min(HH_H_sort)) / (max(HH_H_sort) + min(HH_H_sort));
vis_D = (max(HH_D_sort) - min(HH_D_sort)) / (max(HH_D_sort) + min(HH_D_sort));

fprintf('Alice at H --> Visibility: %.2f%%', vis_H*100);
disp(' ');
fprintf('Alice at D --> Visibility: %.2f%%', vis_D*100);



%% Plotting
figure('Name','Interference Fringes','Color','w','Position',[100 100 1200 500]);

% -------- Plot 1: H --------
ax1 = subplot(1,2,1);
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

if ~isempty(ang_H_sort)
    xq_H = linspace(min(ang_H_sort), max(ang_H_sort), 200);

    % curve (same as your working code):
    yq_H_HH = spline(ang_H_sort, HH_H_sort, xq_H);
    yq_H_HV = spline(ang_H_sort, HV_H_sort, xq_H);

    plot(ax1, xq_H, yq_H_HH, '-', ...
        'Color',[0.35 0.60 0.90], 'LineWidth',1.5, 'HandleVisibility','off');

    plot(ax1, xq_H, yq_H_HV, '-', ...
        'Color',[0.92, 0.45, 0.45], 'LineWidth',1.5, 'HandleVisibility','off');

    plot(ax1, ang_H_sort, HH_H_sort, 'o', ...
        'MarkerFaceColor',[0.35 0.60 0.90], ...
        'MarkerEdgeColor','k', 'MarkerSize',7, 'DisplayName','HH');

    plot(ax1, ang_H_sort, HV_H_sort, 'o', ...
        'MarkerFaceColor',[0.92, 0.45, 0.45], ...
        'MarkerEdgeColor','k', 'MarkerSize',7, 'DisplayName','HV');
end

xlabel(ax1, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize',12);
ylabel(ax1, 'True Coincidences (cps)', 'FontSize',12);
title(ax1, sprintf('Alice at $H$ | $V = %.1f\\%%$', vis_H*100), 'FontSize',13);
legend(ax1,'Location','best','FontSize',10);
set(ax1,'FontSize',11);

%%
% -------- Plot 2: D --------
ax2 = subplot(1,2,2);
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

if ~isempty(ang_D_sort)
    xq_D = linspace(min(ang_D_sort), max(ang_D_sort), 200);

    yq_D_HH = spline(ang_D_sort, HH_D_sort, xq_D);
    yq_D_HV = spline(ang_D_sort, HV_D_sort, xq_D);


    plot(ax2, xq_D, yq_D_HH, '-', ...
        'Color',[0.35 0.60 0.90], 'LineWidth',1.5, 'HandleVisibility','off');
    plot(ax2, xq_D, yq_D_HV, '-', ...
        'Color',[0.92, 0.45, 0.45], 'LineWidth',1.5, 'HandleVisibility','off');

    plot(ax2, ang_D_sort, HH_D_sort, 'o', ...
        'Color',[0.35 0.60 0.90], 'LineWidth',1.5, ...
        'MarkerSize',7, 'DisplayName','HH');
    plot(ax2, ang_D_sort, HV_D_sort, 'o', ...
        'MarkerFaceColor',[0.92, 0.45, 0.45], ...
        'MarkerEdgeColor','k', 'MarkerSize',7, 'DisplayName','HV');
end

xlabel(ax2, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize',12);
ylabel(ax2, 'True Coincidences (cps)', 'FontSize',12);
title(ax2, sprintf('Alice at $D$ | $V = %.1f\\%%$', vis_D*100), 'FontSize',13);
legend(ax2,'Location','best','FontSize',10);
set(ax2,'FontSize',11);

%% sin fit
%{
% Define cosine model: y = a*cos(b*x + c) + d
cosineModel = @(p, x) p(1) * cos(p(2) * x + p(3)) + p(4);
 
% Initial parameter guesses [a, b, c, d]
p0 = [max(y), 1, 0, mean(y)];
 
% Fit the model
pFit = lsqcurvefit(cosineModel, p0, x, y);
 
% Plot results
plot(x, y, 'bo', x, cosineModel(pFit, x), 'r-', 'LineWidth', 2);
legend('Data', 'Cosine Fit');
%}

% x = ang_H_sort; y = HH_H_sort;
% cosineModel = @(p, x) p(1) * cos(p(2) * x + p(3)) + p(4);
% p0 = [max(y), pi/360, 0, mean(y)];
% pFit = lsqcurvefit(cosineModel, p0, x, y);
% plot(x, y, 'bo', x, cosineModel(pFit, x), 'r-', 'LineWidth', 2);
% legend('Data', 'Cosine Fit');