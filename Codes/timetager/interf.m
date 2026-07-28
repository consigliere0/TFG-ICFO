% =========================================================================
%% TFG - ENTANGLEMENT INTERFERENCE FRINGES & VISIBILITY
% =========================================================================
clear; close all; format long g

% -------------------------------------------------------------------------
% 1. INITIAL SETTINGS
% -------------------------------------------------------------------------
set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

data_path = '../../Tests/timetagger/wvg2/entr2'; % Canvia-ho a la ruta de la teva carpeta si cal
tau = 400 * 1e-12; % Finestra de coincidència (400 ps)

% -------------------------------------------------------------------------
% 2. FIND AND PARSE FILES
% -------------------------------------------------------------------------
files = dir(fullfile(data_path, '*INTERFERENCE_A-*_B-*.txt'));

% Arrays per guardar les dades
angles_H = []; trueHH_H = []; trueHV_H = [];
angles_D = []; trueHH_D = []; trueHV_D = [];

disp('Processant arxius...');

for i = 1:length(files)
    fname = files(i).name;
    full_path = fullfile(files(i).folder, fname);
    
    % Extreure la base d'Alice i l'angle de Bob del nom de l'arxiu
    tokens = regexp(fname, 'INTERFERENCE_A-(H|D)_B-(\d{3})', 'tokens');
    if isempty(tokens)
        continue;
    end
    
    alice_base = tokens{1}{1};
    bob_angle  = str2double(tokens{1}{2});
    
    % LLegir dades esquivant l'encapçalament mal format (saltem 12 línies)
    try
        data = readmatrix(full_path, 'NumHeaderLines', 12);
    catch
        warning('No s''ha pogut llegir: %s', fname); continue;
    end
    
    % Si l'arxiu està buit, saltar-lo
    if isempty(data)
        continue;
    end

    % L'arxiu només té 1 fila de dades, així que la mitjana és trivial
    % Índexs basats en la línia: 0, timestamp(text->NaN), CH1, CH2, CH3, CH4, CH5, C1, C2, C3
    S_HA = data(1, 3);
    S_HB = data(1, 5);
    S_VB = data(1, 6);
    
    raw_HH = data(1, 8);  % Coinc1 (HH)
    raw_HV = data(1, 10); % Coinc3 (HV) -> Ara a la columna 10!
    
    % Càlcul de coincidències reals
    acc_HH = S_HA * S_HB * tau;
    acc_HV = S_HA * S_VB * tau;
    
    t_HH = raw_HH - acc_HH;
    t_HV = raw_HV - acc_HV;
    
    % Separar dades segons l'estat d'Alice
    if strcmp(alice_base, 'H')
        angles_H(end+1) = bob_angle;
        trueHH_H(end+1) = t_HH;
        trueHV_H(end+1) = t_HV;
    elseif strcmp(alice_base, 'D')
        angles_D(end+1) = bob_angle;
        trueHH_D(end+1) = t_HH;
        trueHV_D(end+1) = t_HV;
    end
end

% -------------------------------------------------------------------------
% 3. SORT DATA & HANDLE DUPLICATE ANGLES
% -------------------------------------------------------------------------
% H basis: Find unique angles, sort them, and average any duplicate measurements
if ~isempty(angles_H)
    [angles_H, ~, idx_H] = unique(angles_H);
    % L'operador (:) força matrius de columna Nx1 per evitar errors de dimensió
    trueHH_H = accumarray(idx_H(:), trueHH_H(:), [], @mean)';
    trueHV_H = accumarray(idx_H(:), trueHV_H(:), [], @mean)';
end

% D basis: Find unique angles, sort them, and average any duplicate measurements
if ~isempty(angles_D)
    [angles_D, ~, idx_D] = unique(angles_D);
    trueHH_D = accumarray(idx_D(:), trueHH_D(:), [], @mean)';
    trueHV_D = accumarray(idx_D(:), trueHV_D(:), [], @mean)';
end
% -------------------------------------------------------------------------
% 4. CALCULATE VISIBILITY & INFERRED CHSH (S)
% -------------------------------------------------------------------------
if ~isempty(trueHH_H)
    max_H = max(trueHH_H);
    min_H = min(trueHH_H);
    vis_H = (max_H - min_H) / (max_H + min_H);
    %S_inf_H = 2 * sqrt(2) * vis_H;
else
    vis_H = 0; %S_inf_H = 0;
end

if ~isempty(trueHH_D)
    max_D = max(trueHH_D);
    min_D = min(trueHH_D);
    vis_D = (max_D - min_D) / (max_D + min_D);
    %S_inf_D = 2 * sqrt(2) * vis_D;
else
    vis_D = 0; %S_inf_D = 0;
end

fprintf('\n--- RESULTS ---\n');
%fprintf('Alice at H -> Visibility: %.2f%% | Inferred S: %.3f\n', vis_H*100, S_inf_H);
%fprintf('Alice at D -> Visibility: %.2f%% | Inferred S: %.3f\n', vis_D*100, S_inf_D);
fprintf('Alice at H -> Visibility: %.2f%%', vis_H*100);
fprintf('Alice at D -> Visibility: %.2f%%', vis_D*100);

% -------------------------------------------------------------------------
% 5. PLOTTING
% -------------------------------------------------------------------------
if isempty(angles_H) && isempty(angles_D)
    error('No s''han pogut extreure dades dels arxius. Revisa la ruta o els noms dels fitxers.');
end

figure('Name', 'Interference Fringes', 'Color', 'w', 'Position', [100 100 1200 500]);

% --- PLOT 1: Alice at H ---
ax1 = subplot(1, 2, 1);
hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');

if ~isempty(angles_H)
    xq_H = linspace(min(angles_H), max(angles_H), 200);
    yq_HH_H = spline(angles_H, trueHH_H, xq_H);
    yq_HV_H = spline(angles_H, trueHV_H, xq_H);
    
    plot(ax1, xq_H, yq_HH_H, '-', 'Color', [0.35, 0.60, 0.90], 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax1, angles_H, trueHH_H, 'o', 'MarkerFaceColor', [0.35, 0.60, 0.90], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'HH');
    
    plot(ax1, xq_H, yq_HV_H, '-', 'Color', [0.92, 0.45, 0.45], 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax1, angles_H, trueHV_H, 's', 'MarkerFaceColor', [0.92, 0.45, 0.45], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'HV');
end

xlabel(ax1, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax1, 'True Coincidences (cps)', 'FontSize', 12);
%title(ax1, sprintf('Alice at $H | V = %.1f\\%%$ | $S_{inf} = %.2f$', vis_H*100, S_inf_H), 'FontSize', 13);
title(ax1, sprintf('Alice at $H$ | $V = %.1f\\%%$', vis_H*100), 'FontSize', 13);
legend(ax1, 'Location', 'best', 'FontSize', 10);
set(ax1, 'FontSize', 11);

% --- PLOT 2: Alice at D ---
ax2 = subplot(1, 2, 2);
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');

if ~isempty(angles_D)
    xq_D = linspace(min(angles_D), max(angles_D), 200);
    yq_HH_D = spline(angles_D, trueHH_D, xq_D);
    yq_HV_D = spline(angles_D, trueHV_D, xq_D);
    
    plot(ax2, xq_D, yq_HH_D, '-', 'Color', [0.35, 0.60, 0.90], 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax2, angles_D, trueHH_D, 'o', 'MarkerFaceColor', [0.35, 0.60, 0.90], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'HH');
    
    plot(ax2, xq_D, yq_HV_D, '-', 'Color', [0.92, 0.45, 0.45], 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax2, angles_D, trueHV_D, 's', 'MarkerFaceColor', [0.92, 0.45, 0.45], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'HV');
end

xlabel(ax2, 'Bob''s Angle $\theta_{B}$ (deg)', 'FontSize', 12);
ylabel(ax2, 'True Coincidences (cps)', 'FontSize', 12);
%title(ax2, sprintf('Alice at $D$ | $V = %.1f\\%% $| $S_{inf} = %.2f$', vis_D*100, S_inf_D), 'FontSize', 13);
title(ax2, sprintf('Alice at $D$ | $V = %.1f\\%%$', vis_D*100), 'FontSize', 13);
legend(ax2, 'Location', 'best', 'FontSize', 10);
set(ax2, 'FontSize', 11);

%sgtitle('Quantum Interference Fringes', 'FontSize', 15, 'Interpreter', 'latex');