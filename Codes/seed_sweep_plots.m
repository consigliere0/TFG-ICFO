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
dataFolder = '../Tests/Attenuator_Sweep_Test_02';
filePattern = fullfile(dataFolder, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);

% 1. Constants del Sistema
c = 299792458;                   
wl_pump = 1549.65e-9; % La posició del Pump "vista" per l'OSA
wl_seed = 1545.00e-9; % La posició fixa que hem configurat al Seed

% Calcular on està l'Idler (Serà un punt fix per a tot l'experiment)
f_pump = c / wl_pump;
f_seed = c / wl_seed;
f_idler = 2 * f_pump - f_seed;
wl_idler = c / f_idler; 

search_window = 0.5e-9; % Finestra de +/- 0.5 nm

% Pre-assignació de memòria
P_seed_dBm = zeros(numFiles, 1);
P_idler_dBm = zeros(numFiles, 1);

% --- EXTRACCIÓ DE DADES ---
for k = 1:numFiles
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    data = readmatrix(fullFileName);
    
    x = data(:,1); % Longitud d'ona
    y = data(:,2); % Potència en dBm
    
    % 1. Trobar la potència real del Seed (Prop de wl_seed)
    idx_seed = (x >= (wl_seed - search_window)) & (x <= (wl_seed + search_window));
    P_seed_dBm(k) = max(y(idx_seed));
    
    % 2. Trobar la potència real de l'Idler generat (Prop de wl_idler)
    idx_idler = (x >= (wl_idler - search_window)) & (x <= (wl_idler + search_window));
    % Failsafe per si el soroll tapa l'Idler a potències molt baixes
    if any(idx_idler)
        P_idler_dBm(k) = max(y(idx_idler));
    else
        P_idler_dBm(k) = NaN;
    end
end

% Sort function: Com que l'ordre dels fitxers llegits per 'dir' pot no ser seqüencial 
% (ex: llegeix el 10.00dBm abans que el 2.00dBm), ordenem els vectors respecte el Seed.
[P_seed_dBm, sortIdx] = sort(P_seed_dBm);
P_idler_dBm = P_idler_dBm(sortIdx);

% --- CONVERSIÓ A LINEAL (Watts a mW) ---
% Formula: P(mW) = 10^(P(dBm) / 10)
P_seed_mW = 10.^(P_seed_dBm / 10);
P_idler_mW = 10.^(P_idler_dBm / 10);

% --- VISUALITZACIÓ DEL RESULTAT ---
fig = figure('Name', 'Idler Power vs Seed Power', 'Position', [100, 100, 900, 400]);
sgtitle('Idler Power vs Seed Power')

% --- SUBPLOT 1: Escala Lineal (mW) ---
% Hauries de veure una línia recta que passa per l'origen (y = mx)
subplot(1, 2, 1);
plot(P_seed_mW, P_idler_mW, 'squarek', 'LineWidth', 1.5, 'MarkerFaceColor', '#000000');
grid on;
xlabel('Seed Power (mW)');
ylabel('Idler Power (mW)');
title('Linear scale (mW)')

% Opcional: Afegir un fit lineal per comprovar la teoria
p = polyfit(P_seed_mW, P_idler_mW, 1);
yfit = polyval(p, P_seed_mW);
hold on; plot(P_seed_mW, yfit, '-', 'LineWidth', 1, 'Color', '#a9a9a9'); hold off;


% --- SUBPLOT 2: Escala Logarítmica (dBm) ---
% Hauries de veure una línia recta amb pendent = 1
subplot(1, 2, 2);
plot(P_seed_dBm, P_idler_dBm, 'pentagramk', 'LineWidth', 1.5, 'MarkerFaceColor', '#000000');
grid on;
xlabel('Seed Power (dBm)');
ylabel('Idler Power (dBm)');
title('Logarithmic scale (dBm)');

% Opcional: Afegir un fit lineal per comprovar la teoria
p = polyfit(P_seed_dBm, P_idler_dBm, 1);
yfit = polyval(p, P_seed_dBm);
hold on; plot(P_seed_dBm, yfit, '-', 'LineWidth', 1, 'Color', '#a9a9a9'); 

% Add slope label
slope_val = p(1);
txt = sprintf('m = %.3g dB/dB', slope_val);
xpos = -34; %mean(P_seed_dBm(valid));
ypos = -51;mean(P_idler_dBm);
t = text(xpos, ypos, txt,'BackgroundColor','white','EdgeColor','k');
uistack(t,'top');
hold off

outFolder = '../Tests/Attenuator_Sweep_Test_02';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
saveas(fig, fullfile(outFolder, 'Seed_power_sweep.png'));





% --- SUBPLOT 2: Escala Logarítmica (dBm) ---
% Hauries de veure una línia recta amb pendent = 1
% subplot(1,2,2);
% valid = ~isnan(P_idler_dBm) & ~isnan(P_seed_dBm);
% plot(P_seed_dBm(valid), P_idler_dBm(valid), 'pentagramk', 'LineWidth', 1.5, 'MarkerFaceColor', '#000000');
% grid on;
% xlabel('Seed Power (dBm)');
% ylabel('Idler Power (dBm)');
% title('Logarithmic scale (dBm)');
% 
% p = polyfit(P_seed_dBm(valid), P_idler_dBm(valid), 1);
% yfit = polyval(p, P_seed_dBm(valid));
% hold on;
% plot(P_seed_dBm(valid), yfit, '-', 'LineWidth', 1, 'Color', '#a9a9a9');
% 
% slope_val = p(1);
% txt = sprintf('m = %.3g dB/dB', slope_val);
% xpos = -31; %mean(P_seed_dBm(valid));
% ypos = -53;mean(P_idler_dBm(valid));
% t = text(xpos, ypos, txt,'BackgroundColor','white','EdgeColor','k');
% uistack(t,'top');
% hold off