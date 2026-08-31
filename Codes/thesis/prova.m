% =========================================================================
%%                     ENTANGLEMENT INTERFERENCE + CHSH
% =========================================================================
clear all; close all; format long g; clc;

set(groot, 'defaultTextInterpreter',         'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

% *******************************
%%       Data processing
% *******************************
tau = 400 * 1e-12; % coincidence window 400ps
acq_time = 10;     % <--- SET YOUR ACQUISITION TIME PER POINT (SECONDS) HERE

data_path = '../../Tests/timetagger/wvg2/entr9';
files = dir(fullfile(data_path, '*A*-B*.txt')); 

% ARRAYS to hold data
% Procedure 1: HWP at 0 deg -> Alice projects to 0 deg (Ch1) and 90 deg (Ch2)
ang_B_1 = [];
A0_B    = []; A0_B90  = []; % Alice Transmitted (0 deg) coincidences
A90_B   = []; A90_B90 = []; % Alice Reflected (90 deg) coincidences

% Procedure 2: HWP at 22.5 deg -> Alice projects to 45 deg (Ch1) and -45 deg (Ch2)
ang_B_2 = [];
A45_B   = []; A45_B90  = []; % Alice Transmitted (45 deg) coincidences
Am45_B  = []; Am45_B90 = []; % Alice Reflected (-45 deg) coincidences

for k = 1:length(files)
    file_path = fullfile(data_path, files(k).name);
    data = readmatrix(file_path, 'NumHeaderLines', 12);  
    
    RM_Alice = data(1, 13); 
    RM_Bob   = data(1, 12) * 2; % Bob's angle delta_B
    
    % Single counts
    ch1 = data(1, 3); 
    ch2 = data(1, 4); 
    ch3 = data(1, 5); 
    ch4 = data(1, 6); 
    
    % Raw Coincidences
    raw_13 = data(1, 8);  % Alice Transmitted & Bob Transmitted
    raw_24 = data(1, 9);  % Alice Reflected & Bob Reflected
    raw_14 = data(1, 10); % Alice Transmitted & Bob Reflected
    raw_23 = data(1, 11); % Alice Reflected & Bob Transmitted
    
    % Accidentals calculation
    acc_13 = ch1 * ch3 * tau;
    acc_14 = ch1 * ch4 * tau;
    acc_24 = ch2 * ch4 * tau;
    acc_23 = ch2 * ch3 * tau;
    
    if RM_Alice == 0.0 % Procedure 1
        ang_B_1 = [ang_B_1, RM_Bob];
        A0_B    = [A0_B,    raw_13 - acc_13];
        A0_B90  = [A0_B90,  raw_14 - acc_14];
        A90_B   = [A90_B,   raw_23 - acc_23];
        A90_B90 = [A90_B90, raw_24 - acc_24];
    elseif RM_Alice == 22.5 % Procedure 2
        ang_B_2 = [ang_B_2, RM_Bob];
        A45_B    = [A45_B,    raw_13 - acc_13];
        A45_B90  = [A45_B90,  raw_14 - acc_14];
        Am45_B   = [Am45_B,   raw_23 - acc_23];
        Am45_B90 = [Am45_B90, raw_24 - acc_24];
    end
end

% Sort Procedure 1
[ang_B1_sort, idx_1] = sort(ang_B_1, 'ascend');
A0_B   = A0_B(idx_1);
A0_B90 = A0_B90(idx_1);
A90_B  = A90_B(idx_1);
A90_B90= A90_B90(idx_1);

% Sort Procedure 2
[ang_B2_sort, idx_2] = sort(ang_B_2, 'ascend');
A45_B   = A45_B(idx_2);
A45_B90 = A45_B90(idx_2);
Am45_B  = Am45_B(idx_2);
Am45_B90= Am45_B90(idx_2);

% Poisson Error Calculation: sqrt(True_Counts) / T 
err_A0_B   = sqrt(max(A0_B, 0)   * acq_time) / acq_time;
err_A0_B90 = sqrt(max(A0_B90, 0) * acq_time) / acq_time;
err_A90_B  = sqrt(max(A90_B, 0)  * acq_time) / acq_time;
err_A90_B90= sqrt(max(A90_B90, 0)* acq_time) / acq_time;

err_A45_B   = sqrt(max(A45_B, 0)   * acq_time) / acq_time;
err_A45_B90 = sqrt(max(A45_B90, 0) * acq_time) / acq_time;
err_Am45_B  = sqrt(max(Am45_B, 0)  * acq_time) / acq_time;
err_Am45_B90= sqrt(max(Am45_B90, 0)* acq_time) / acq_time;

% *******************************
%%       Cosine fit & Visibility
% *******************************
cosineModel = @(p, x) p(1) * cos(p(2) * x + p(3)) + p(4);
opts = optimoptions('lsqcurvefit', 'Display', 'off');
fitCosine = @(x, y) local_fit_cosine(x, y, cosineModel, opts);

% Fit Procedure 1
[p_A0_B,  xq_1, yq_A0_B]   = fitCosine(ang_B1_sort, A0_B);
[p_A0_B90,~,    yq_A0_B90] = fitCosine(ang_B1_sort, A0_B90);
[p_A90_B, ~,    yq_A90_B]  = fitCosine(ang_B1_sort, A90_B);
[p_A90_B90,~,   yq_A90_B90]= fitCosine(ang_B1_sort, A90_B90);

% Fit Procedure 2
[p_A45_B,  xq_2, yq_A45_B]   = fitCosine(ang_B2_sort, A45_B);
[p_A45_B90,~,    yq_A45_B90] = fitCosine(ang_B2_sort, A45_B90);
[p_Am45_B, ~,    yq_Am45_B]  = fitCosine(ang_B2_sort, Am45_B);
[p_Am45_B90,~,   yq_Am45_B90]= fitCosine(ang_B2_sort, Am45_B90);

% Visibilities: V = |a| / d
v_A0_B   = abs(p_A0_B(1))   / p_A0_B(4);
v_A0_B90 = abs(p_A0_B90(1)) / p_A0_B90(4);
v_A90_B  = abs(p_A90_B(1))  / p_A90_B(4);
v_A90_B90= abs(p_A90_B90(1))/ p_A90_B90(4);

v_A45_B   = abs(p_A45_B(1))   / p_A45_B(4);
v_A45_B90 = abs(p_A45_B90(1)) / p_A45_B90(4);
v_Am45_B  = abs(p_Am45_B(1))  / p_Am45_B(4);
v_Am45_B90= abs(p_Am45_B90(1))/ p_Am45_B90(4);

overall_vis_0_90 = mean([v_A0_B, v_A0_B90, v_A90_B, v_A90_B90]);
overall_vis_45_m45 = mean([v_A45_B, v_A45_B90, v_Am45_B, v_Am45_B90]);

fprintf('--- INTERFERENCE VISIBILITIES ---\n');
fprintf('Overall Rectilinear Basis Visibility: %.2f%%\n', overall_vis_0_90 * 100);
fprintf('Overall Diagonal Basis Visibility: %.2f%%\n', overall_vis_45_m45 * 100);

% ********************************************
%%       Plotting 
% ********************************************
figure('Name','Interference Fringes Phase Shift','Color','w','Position',[100 100 900 700]);

% Tiled layout for shared global axes 
t = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

%col_R = [0.00, 0.45, 0.74]; 
%col_D = [0.85, 0.33, 0.10]; 

col_R =    '#4CAF50'; % Loser TM: Verd pastel fosc
col_D =  '#7E57C2'; % TM: lila fosc

% --- Subplot 1: Ch1 & Ch3 (Alice 1 & Bob 1) ---
ax1 = nexttile(1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, xq_1, yq_A0_B,  '-', 'Color', col_R, 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax1, xq_2, yq_A45_B, '-', 'Color', col_D, 'LineWidth', 1.5, 'HandleVisibility','off');

errorbar(ax1, ang_B1_sort, A0_B, err_A0_B, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_R, ...
    'Color', col_R, 'DisplayName', '$\delta_A = 0^\circ$');
errorbar(ax1, ang_B2_sort, A45_B, err_A45_B, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_D, ...
    'Color', col_D, 'DisplayName', '$\delta_A = 45^\circ$');

% Text box mimicking the paper's style
text(ax1, 0.95, 0.95, {'Alice 1', 'Bob 1'}, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'k', 'BackgroundColor', [0.9 0.9 0.9], 'Margin', 3, 'Interpreter', 'latex');

%xticklabels(ax1, {}); % Hide x-axis tick labels for top plots
xlim(ax1, [0,360]); set(ax1, 'FontSize', 11);

% --- Subplot 2: Ch1 & Ch4 (Alice 1 & Bob 2) ---
ax2 = nexttile(2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, xq_1, yq_A0_B90,  '-', 'Color', col_R, 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax2, xq_2, yq_A45_B90, '-', 'Color', col_D, 'LineWidth', 1.5, 'HandleVisibility','off');

errorbar(ax2, ang_B1_sort, A0_B90, err_A0_B90, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_R, ...
    'Color', col_R, 'HandleVisibility', 'off');
errorbar(ax2, ang_B2_sort, A45_B90, err_A45_B90, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_D, ...
    'Color', col_D, 'HandleVisibility', 'off');

text(ax2, 0.95, 0.95, {'Alice 1', 'Bob 2'}, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'k', 'BackgroundColor', [0.9 0.9 0.9], 'Margin', 3, 'Interpreter', 'latex');

%xticklabels(ax2, {}); % Hide x-axis tick labels
%yticklabels(ax2, {}); % Hide y-axis tick labels for right plots
xlim(ax2, [0,360]); set(ax2, 'FontSize', 11);

% --- Subplot 3: Ch2 & Ch3 (Alice 2 & Bob 1) ---
ax3 = nexttile(3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
plot(ax3, xq_1, yq_A90_B,  '-', 'Color', col_R, 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax3, xq_2, yq_Am45_B, '-', 'Color', col_D, 'LineWidth', 1.5, 'HandleVisibility','off');

errorbar(ax3, ang_B1_sort, A90_B, err_A90_B, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_R, ...
    'Color', col_R, 'HandleVisibility', 'off');
errorbar(ax3, ang_B2_sort, Am45_B, err_Am45_B, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_D, ...
    'Color', col_D, 'HandleVisibility', 'off');

text(ax3, 0.95, 0.95, {'Alice 2', 'Bob 1'}, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'k', 'BackgroundColor', [0.9 0.9 0.9], 'Margin', 3, 'Interpreter', 'latex');

xlim(ax3, [0,360]); set(ax3, 'FontSize', 11);

% --- Subplot 4: Ch2 & Ch4 (Alice 2 & Bob 2) ---
ax4 = nexttile(4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
plot(ax4, xq_1, yq_A90_B90,  '-', 'Color', col_R, 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax4, xq_2, yq_Am45_B90, '-', 'Color', col_D, 'LineWidth', 1.5, 'HandleVisibility','off');

errorbar(ax4, ang_B1_sort, A90_B90, err_A90_B90, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_R, ...
    'Color', col_R, 'HandleVisibility', 'off');
errorbar(ax4, ang_B2_sort, Am45_B90, err_Am45_B90, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col_D, ...
    'Color', col_D, 'HandleVisibility', 'off');

text(ax4, 0.95, 0.95, {'Alice 2', 'Bob 2'}, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'k', 'BackgroundColor', [0.9 0.9 0.9], 'Margin', 3, 'Interpreter', 'latex');

%yticklabels(ax4, {}); % Hide y-axis tick labels
xlim(ax4, [0,360]); set(ax4, 'FontSize', 11);

% Global Axes Labels
xlabel(t, 'Bob''s Angle $\delta_{B}$ (deg)', 'FontSize', 13, 'Interpreter', 'latex');
ylabel(t, 'True Coincidences (cps)', 'FontSize', 13, 'Interpreter', 'latex');

% Single Legend at the top (Matches Figure 5 layout)
lgd = legend(ax1, 'Orientation', 'horizontal', 'FontSize', 12, 'Interpreter', 'latex');
lgd.Layout.Tile = 'North'; 

% *******************************
%%       CHSH: three angle choices, side by side
% *******************************

% Normalized joint correlators E(theta_A, theta_B)
% Procedure 1 -> Alice fixed at 0 deg;  Procedure 2 -> Alice fixed at 45 deg
E1_raw = (A0_B + A90_B90 - A0_B90 - A90_B) ./ (A0_B + A90_B90 + A0_B90 + A90_B);      % E(0,  theta_B)
E2_raw = (A45_B + Am45_B90 - A45_B90 - Am45_B) ./ (A45_B + Am45_B90 + A45_B90 + Am45_B); % E(45, theta_B)

% Fit each correlator vs Bob's angle so we can evaluate E at ANY angle,
% not just the ones actually measured
[pE1, ~, ~] = fitCosine(ang_B1_sort, E1_raw);
[pE2, ~, ~] = fitCosine(ang_B2_sort, E2_raw);
E1fun = @(th) cosineModel(pE1, th);   % E(0,  theta_B), continuous
E2fun = @(th) cosineModel(pE2, th);   % E(45, theta_B), continuous

% % ---- (1) Literal supervisor angles: Bob at 0 / 45 ----
% [S1, sS1] = chsh_at(0, 45, ang_B1_sort, ang_B2_sort, A0_B,A0_B90,A90_B,A90_B90, ...
%     A45_B,A45_B90,Am45_B,Am45_B90, err_A0_B,err_A0_B90,err_A90_B,err_A90_B90, ...
%     err_A45_B,err_A45_B90,err_Am45_B,err_Am45_B90, E1fun, E2fun);

% ---- (2) Textbook optimal (assumes exactly 45 deg phase mismatch): Bob at 22.5 / 67.5 ----
[S2, sS2] = chsh_at(22.5, 67.5, ang_B1_sort, ang_B2_sort, A0_B,A0_B90,A90_B,A90_B90, ...
    A45_B,A45_B90,Am45_B,Am45_B90, err_A0_B,err_A0_B90,err_A90_B,err_A90_B90, ...
    err_A45_B,err_A45_B90,err_Am45_B,err_Am45_B90, E1fun, E2fun);

% ---- (3) True optimal from grid search over the fitted curves ----
% (accounts for the real ~55 deg phase mismatch measured on this hardware)
bgrid = 0:0.05:360;
[B, Bp] = meshgrid(bgrid, bgrid);
% NOTE: meshgrid(bgrid,bgrid) gives B(i,j)=bgrid(j), Bp(i,j)=bgrid(i) --
% i.e. the FIRST index (i) corresponds to "bp" and the SECOND (j) to "b".
% Getting this backwards silently produces a different (wrong) S, since
% S is not symmetric under swapping b <-> bp.
Sgrid = E1fun(B) - E1fun(Bp) + E2fun(B) + E2fun(Bp);
[Sgrid_best, idxMax] = max(abs(Sgrid(:)));
[i, j] = ind2sub(size(Sgrid), idxMax);
b_opt  = bgrid(j);   % B(i,j)  = bgrid(j)
bp_opt = bgrid(i);   % Bp(i,j) = bgrid(i)

[S3, sS3] = chsh_at(b_opt, bp_opt, ang_B1_sort, ang_B2_sort, A0_B,A0_B90,A90_B,A90_B90, ...
    A45_B,A45_B90,Am45_B,Am45_B90, err_A0_B,err_A0_B90,err_A90_B,err_A90_B90, ...
    err_A45_B,err_A45_B90,err_Am45_B,err_Am45_B90, E1fun, E2fun);

% Self-check: chsh_at's recomputed S must match the grid value at (i,j),
% otherwise the index mapping above is wrong again
assert(abs(abs(S3) - Sgrid_best) < 1e-9, ...
    'Grid-search index mapping mismatch: S3 does not match Sgrid(idxMax). Check meshgrid orientation.');

% --- Also compute S directly from the raw counts at the NEAREST measured
% Bob angle (supervisor's "method 1"), as opposed to S1/S2/S3 above which
% use the smooth cosine fit ("method 2"). ---
S_raw_at = @(b, bp) local_S_from_raw(b, bp, ang_B1_sort, ang_B2_sort, ...
    A0_B,A0_B90,A90_B,A90_B90, A45_B,A45_B90,Am45_B,Am45_B90);

S1_raw = S_raw_at(0, 45);
S2_raw = S_raw_at(22.5, 67.5);
S3_raw = S_raw_at(b_opt, bp_opt);

fprintf('\n--- CHSH: comparison of angle choices ---\n');
fprintf('%-32s %10s %10s %10s %10s %10s\n', 'Setting', 'S_fit', 'S_raw', 'sigma_S', 'N_sigma(fit)', 'N_sigma(raw)');
%fprintf('%-32s %10.4f %10.4f %10.4f %10.2f %10.2f\n', 'Literal 0/45', S1, S1_raw, sS1, (abs(S1)-2)/sS1, (abs(S1_raw)-2)/sS1);
fprintf('%-32s %10.4f %10.4f %10.4f %10.2f %10.2f\n', 'Textbook 22.5/67.5', S2, S2_raw, sS2, (abs(S2)-2)/sS2, (abs(S2_raw)-2)/sS2);
fprintf('%-32s %10.4f %10.4f %10.4f %10.2f %10.2f\n', sprintf('True optimal (%.1f/%.1f)', b_opt, bp_opt), S3, S3_raw, sS3, (abs(S3)-2)/sS3, (abs(S3_raw)-2)/sS3);
fprintf('\n(S_fit: correlators evaluated on the smooth cosine fit. S_raw: correlators\n');
fprintf(' taken from the raw coincidence counts nearest each target angle.\n');
fprintf(' sigma_S is the Poisson-propagated uncertainty from the raw counts, used for both.)\n');


% =========================================================================
%%          Exact phase offset from the fits (no eyeballing) 
% =========================================================================
w1 = pE1(2); phi1 = pE1(3);
w2 = pE2(2); phi2 = pE2(3);
period1 = 2*pi/w1; period2 = 2*pi/w2;   % should both be ~180 deg

theta_peak_1 = mod(-phi1/w1, period1);   % Alice=H: expect peak near Bob=0
theta_peak_2 = mod(-phi2/w2, period2);   % Alice=D: expect peak near Bob=45

raw_offset = theta_peak_2 - theta_peak_1;
measured_offset_deg = raw_offset - 180*round((raw_offset - 45)/180);  % wrap near 45

fprintf('\n--- Measured phase offset (exact) ---\n');
fprintf('Peak of E(0,theta_B)  at %.2f deg (period %.2f deg)\n', theta_peak_1, period1);
fprintf('Peak of E(45,theta_B) at %.2f deg (period %.2f deg)\n', theta_peak_2, period2);
fprintf('Measured offset: %.2f deg (theory expects exactly 45 deg)\n', measured_offset_deg);
fprintf('Excess beyond theory: %.2f deg\n', measured_offset_deg - 45);

% Independent cross-check using the RAW fringe fits (not the E-correlator
% fits), for consistency
w1b=p_A0_B(2); phi1b=p_A0_B(3); w2b=p_A45_B(2); phi2b=p_A45_B(3);
theta_peak_1b = mod(-phi1b/w1b, 2*pi/w1b);
theta_peak_2b = mod(-phi2b/w2b, 2*pi/w2b);
offset_b = theta_peak_2b - theta_peak_1b;
measured_offset_deg_b = offset_b - 180*round((offset_b - 45)/180);
fprintf('Cross-check from raw A0_B/A45_B fits: %.2f deg\n', measured_offset_deg_b);

%% --- Constrained bisector search: b, b' forced exactly 45 deg apart ---
betaGrid = 0:0.02:180;
S_of_beta = E1fun(betaGrid-22.5) - E1fun(betaGrid+22.5) ...
          + E2fun(betaGrid-22.5) + E2fun(betaGrid+22.5);
[~, iBeta] = max(abs(S_of_beta));
beta_opt = betaGrid(iBeta);
b_con  = beta_opt - 22.5;
bp_con = beta_opt + 22.5;

fprintf('\n--- Constrained bisector search ---\n');
fprintf('Optimal bisector = %.2f deg -> b=%.2f, bp=%.2f\n', beta_opt, b_con, bp_con);
fprintf('Theory-predicted bisector (22.5 + offset/2) = %.2f deg\n', 22.5 + measured_offset_deg/2);

[S_con, sS_con] = chsh_at(b_con, bp_con, ang_B1_sort, ang_B2_sort, A0_B,A0_B90,A90_B,A90_B90, ...
    A45_B,A45_B90,Am45_B,Am45_B90, err_A0_B,err_A0_B90,err_A90_B,err_A90_B90, ...
    err_A45_B,err_A45_B90,err_Am45_B,err_Am45_B90, E1fun, E2fun);
S_con_raw = S_raw_at(b_con, bp_con);

fprintf('S (fit)=%.4f  S (raw)=%.4f  sigma_S=%.4f  N_sigma(raw)=%.2f\n', ...
        S_con, S_con_raw, sS_con, (abs(S_con_raw)-2)/sS_con);



% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================
function [pFit, xq, yq] = local_fit_cosine(x, y, cosineModel, opts)
    x = x(:); y = y(:);
    d0 = mean(y);
    a0 = (max(y) - min(y)) / 2;
    b0 = 2*pi/180;                  
    [~, idx_max] = max(y);
    c0 = -b0 * x(idx_max);
    p0 = [a0, b0, c0, d0];
    
    lb = [-Inf, 0.5*b0, -Inf, -Inf];
    ub = [ Inf, 1.5*b0,  Inf,  Inf];
    pFit = lsqcurvefit(cosineModel, p0, x, y, lb, ub, opts);
    
    xq = linspace(min(x), max(x), 300);
    yq = cosineModel(pFit, xq);
end

function [S, sigma_S, ib, ibp] = chsh_at(b, bp, ang1, ang2, A0_B,A0_B90,A90_B,A90_B90, ...
        A45_B,A45_B90,Am45_B,Am45_B90, err_A0_B,err_A0_B90,err_A90_B,err_A90_B90, ...
        err_A45_B,err_A45_B90,err_Am45_B,err_Am45_B90, E1fun, E2fun)
    % Computes CHSH S at settings (b, bp) using the fitted E-functions,
    % and propagates Poissonian error using the raw coincidence counts
    % nearest to those angles.

    [~, ib1]  = min(abs(ang1 - b));
    [~, ib1p] = min(abs(ang1 - bp));
    [~, ib2]  = min(abs(ang2 - b));
    [~, ib2p] = min(abs(ang2 - bp));

    S = E1fun(b) - E1fun(bp) + E2fun(b) + E2fun(bp);

    sig1  = corr_error(A0_B(ib1),A0_B90(ib1),A90_B(ib1),A90_B90(ib1), ...
                        err_A0_B(ib1),err_A0_B90(ib1),err_A90_B(ib1),err_A90_B90(ib1));
    sig1p = corr_error(A0_B(ib1p),A0_B90(ib1p),A90_B(ib1p),A90_B90(ib1p), ...
                        err_A0_B(ib1p),err_A0_B90(ib1p),err_A90_B(ib1p),err_A90_B90(ib1p));
    sig2  = corr_error(A45_B(ib2),A45_B90(ib2),Am45_B(ib2),Am45_B90(ib2), ...
                        err_A45_B(ib2),err_A45_B90(ib2),err_Am45_B(ib2),err_Am45_B90(ib2));
    sig2p = corr_error(A45_B(ib2p),A45_B90(ib2p),Am45_B(ib2p),Am45_B90(ib2p), ...
                        err_A45_B(ib2p),err_A45_B90(ib2p),err_Am45_B(ib2p),err_Am45_B90(ib2p));

    % Combined in quadrature: valid because these four correlators are
    % measured at four different (independent) Bob settings / photon samples
    sigma_S = sqrt(sig1^2 + sig1p^2 + sig2^2 + sig2p^2);
    ib = b; ibp = bp;
end

function S = local_S_from_raw(b, bp, ang1, ang2, A0_B,A0_B90,A90_B,A90_B90, A45_B,A45_B90,Am45_B,Am45_B90)
    % CHSH S computed directly from the raw coincidence counts nearest to
    % the requested Bob angles b, bp (no fitting involved).
    [~, ib1]  = min(abs(ang1 - b));
    [~, ib1p] = min(abs(ang1 - bp));
    [~, ib2]  = min(abs(ang2 - b));
    [~, ib2p] = min(abs(ang2 - bp));

    E_a0_b   = (A0_B(ib1)  + A90_B90(ib1)  - A0_B90(ib1)  - A90_B(ib1))  / (A0_B(ib1)  + A90_B90(ib1)  + A0_B90(ib1)  + A90_B(ib1));
    E_a0_bp  = (A0_B(ib1p) + A90_B90(ib1p) - A0_B90(ib1p) - A90_B(ib1p)) / (A0_B(ib1p) + A90_B90(ib1p) + A0_B90(ib1p) + A90_B(ib1p));
    E_a45_b  = (A45_B(ib2)  + Am45_B90(ib2)  - A45_B90(ib2)  - Am45_B(ib2))  / (A45_B(ib2)  + Am45_B90(ib2)  + A45_B90(ib2)  + Am45_B(ib2));
    E_a45_bp = (A45_B(ib2p) + Am45_B90(ib2p) - A45_B90(ib2p) - Am45_B(ib2p)) / (A45_B(ib2p) + Am45_B90(ib2p) + A45_B90(ib2p) + Am45_B(ib2p));

    S = E_a0_b - E_a0_bp + E_a45_b + E_a45_bp;
end

function sigE = corr_error(a,b,c,d, sa,sb,sc,sd)
    % Poissonian error propagation for E = (a+d-b-c)/(a+b+c+d)
    % a,b,c,d are independent coincidence counts with std devs sa,sb,sc,sd
    N = a+b+c+d;
    dEda = 2*(b+c)/N^2;   % = dE/dd  (a and d enter symmetrically)
    dEdb = -2*(a+d)/N^2;  % = dE/dc  (b and c enter symmetrically)
    sigE = sqrt( dEda^2*(sa^2+sd^2) + dEdb^2*(sb^2+sc^2) );
end


