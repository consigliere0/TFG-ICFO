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
dataFolder = '../Tests/Pump_Sweep_01'; 
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);

% 1. Constants del Sistema
c = 299792458;                   
wl_pump = 1549.65e-9; % La posició del Pump "vista" per l'OSA
wl_seed = 1545.00e-9; % La posició fixa que hem configurat al Seed

% Calcular on està l'Idler 
f_pump = c / wl_pump;
f_seed = c / wl_seed;
f_idler = 2 * f_pump - f_seed;
wl_idler = c / f_idler; 

search_window = 0.5e-9; % Finestra de +/- 0.5 nm

% Pre-assignació de memòria
P_pump_dBm = zeros(numFiles, 1);
P_idler_dBm = zeros(numFiles, 1);

% --- EXTRACCIÓ DE DADES ---
for k = 1:numFiles
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    data = readmatrix(fullFileName);
    
    x = data(:,1); % wvl
    y = data(:,2); % power dBm
    
    % 1. Find pump power
    idx_pump = (x >= (wl_pump - search_window)) & (x <= (wl_pump + search_window));
    P_pump_dBm(k) = max(y(idx_pump));
    
    % 2. Find exp power of generated idle
    idx_idler = (x >= (wl_idler - search_window)) & (x <= (wl_idler + search_window));
    
    % Failsafe per si el soroll tapa l'Idler a potències molt baixes
    if any(idx_idler)
        P_idler_dBm(k) = max(y(idx_idler));
    else
        P_idler_dBm(k) = NaN;
    end
end

% Sort function: Ordenem els vectors respecte el Pump.
[P_pump_dBm, sortIdx] = sort(P_pump_dBm);
P_idler_dBm = P_idler_dBm(sortIdx);

% Creem màscara de valors vàlids per no fer col·lapsar els 'polyfit'
valid = ~isnan(P_pump_dBm) & ~isnan(P_idler_dBm);
P_pump_dBm_val = P_pump_dBm(valid);
P_idler_dBm_val = P_idler_dBm(valid);

% --- CONVERSIÓ A LINEAL ---
P_pump_mW_val = 10.^(P_pump_dBm_val / 10);
P_idler_mW_val = 10.^(P_idler_dBm_val / 10);

% --- VISUALITZACIÓ  ---
fig = figure('Name', 'Idler Power vs Pump Power', 'Position', [100, 100, 900, 400]);
sgtitle('Idler Power vs Pump Power')

%% PLOT
% --- SUBPLOT 1: Escala Lineal (mW) ---
% Hauries de veure una paràbola (P_idler depèn del quadrat de P_pump)
subplot(1, 2, 1);
plot(P_pump_mW_val, P_idler_mW_val, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', '#000000');
grid on;
xlabel('Pump Power (mW)');
ylabel('Idler Power (mW)');
title('Linear scale (mW)')

% Fit Quadràtic (Grau 2)
p_lin = polyfit(P_pump_mW_val, P_idler_mW_val, 2);
% Generem una x fina per a que la corba es vegi suau
x_fit_lin = linspace(min(P_pump_mW_val), max(P_pump_mW_val), 100);
yfit_lin = polyval(p_lin, x_fit_lin);
hold on; 
plot(x_fit_lin, yfit_lin, '-', 'LineWidth', 1, 'Color', '#a9a9a9'); 
hold off;

% --- SUBPLOT 2: Escala Logarítmica (dBm) ---
% Hauries de veure una línia recta amb pendent ≈ 2
subplot(1, 2, 2);
plot(P_pump_dBm_val, P_idler_dBm_val, 'pentagramk', 'LineWidth', 1.5, 'MarkerFaceColor', '#000000');
grid on;
xlabel('Pump Power (dBm)');
ylabel('Idler Power (dBm)');
title('Logarithmic scale (dBm)');

% Fit Lineal (Grau 1)
p_log = polyfit(P_pump_dBm_val, P_idler_dBm_val, 1);
yfit_log = polyval(p_log, P_pump_dBm_val);
hold on; 
plot(P_pump_dBm_val, yfit_log, '-', 'LineWidth', 1, 'Color', '#a9a9a9'); 

% Add slope label
slope_val = p_log(1);
txt = sprintf('m = %.3g dB/dB', slope_val);
% Ajust automàtic de la posició del text per evitar hardcoding que el deixi fora del gràfic
xpos = min(P_pump_dBm_val) + (max(P_pump_dBm_val)-min(P_pump_dBm_val))*0.1;
ypos = max(P_idler_dBm_val) - (max(P_idler_dBm_val)-min(P_idler_dBm_val))*0.1;
t = text(xpos, ypos, txt, 'BackgroundColor', 'white', 'EdgeColor', 'k');
uistack(t, 'top');
hold off;

% Guardar la figura
outFolder = dataFolder;
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig, fullfile(outFolder, 'Pump_power_sweep.png'));