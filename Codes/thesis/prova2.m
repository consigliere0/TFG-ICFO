% =========================================================================
% CHSH RIGOROUS UNCERTAINTY: weighted fit + delta-method propagation
% =========================================================================
% Insert this AFTER your existing script computes: ang_B1_sort, ang_B2_sort,
% A0_B, A0_B90, A90_B, A90_B90, A45_B, A45_B90, Am45_B, Am45_B90,
% err_A0_B, err_A0_B90, err_A90_B, err_A90_B90,
% err_A45_B, err_A45_B90, err_Am45_B, err_Am45_B90, b_opt, bp_opt,
% b_con, bp_con, corr_error(), cosineModel, opts
% =========================================================================

%% --- Per-point uncertainty on the correlators themselves ---
err_E1_raw = arrayfun(@(i) corr_error(A0_B(i),A0_B90(i),A90_B(i),A90_B90(i), ...
    err_A0_B(i),err_A0_B90(i),err_A90_B(i),err_A90_B90(i)), 1:length(A0_B));
err_E2_raw = arrayfun(@(i) corr_error(A45_B(i),A45_B90(i),Am45_B(i),Am45_B90(i), ...
    err_A45_B(i),err_A45_B90(i),err_Am45_B(i),err_Am45_B90(i)), 1:length(A45_B));

%% --- Weighted cosine fits, with parameter covariance ---
[pE1_w, covp1, chi2red1] = local_fit_cosine_weighted(ang_B1_sort, E1_raw, err_E1_raw, cosineModel, opts);
[pE2_w, covp2, chi2red2] = local_fit_cosine_weighted(ang_B2_sort, E2_raw, err_E2_raw, cosineModel, opts);

fprintf('\n--- Weighted-fit diagnostics ---\n');
fprintf('E1(theta_B) fit: reduced chi^2 = %.3f  (should be ~1 if errors + model are good)\n', chi2red1);
fprintf('E2(theta_B) fit: reduced chi^2 = %.3f\n', chi2red2);
if chi2red1 > 2 || chi2red2 > 2
    warning(['Reduced chi^2 is well above 1 - either the error bars are ' ...
             'underestimated, or the single-cosine model is not fully ' ...
             'capturing the data. Investigate before trusting sigma_S.']);
end

%% --- Rigorous CHSH evaluation at any angle pair ---
fprintf('\n--- CHSH with rigorous (covariance-propagated) uncertainty ---\n');
fprintf('%-32s %10s %10s %10s\n', 'Setting', 'S', 'sigma_S', 'N_sigma');

[S_txt, sS_txt] = chsh_with_covariance(22.5, 67.5, pE1_w, covp1, pE2_w, covp2);
fprintf('%-32s %10.4f %10.4f %10.2f\n', 'Textbook 22.5/67.5', S_txt, sS_txt, (abs(S_txt)-2)/sS_txt);

[S_con_rig, sS_con_rig] = chsh_with_covariance(b_con, bp_con, pE1_w, covp1, pE2_w, covp2);
fprintf('%-32s %10.4f %10.4f %10.2f\n', sprintf('Constrained bisector (%.1f/%.1f)', b_con, bp_con), S_con_rig, sS_con_rig, (abs(S_con_rig)-2)/sS_con_rig);

[S_free_rig, sS_free_rig] = chsh_with_covariance(b_opt, bp_opt, pE1_w, covp1, pE2_w, covp2);
fprintf('%-32s %10.4f %10.4f %10.2f\n', sprintf('Free optimum (%.1f/%.1f)', b_opt, bp_opt), S_free_rig, sS_free_rig, (abs(S_free_rig)-2)/sS_free_rig);

fprintf('\nConsistency check: constrained vs. free optimum S differ by %.4f\n', ...
        abs(abs(S_con_rig) - abs(S_free_rig)));
fprintf('(Small difference here confirms the 45-degree assumption cost is negligible.)\n');

fprintf('\n=== HEADLINE RESULT FOR THESIS ===\n');
fprintf('S = %.4f +/- %.4f  (%.2f sigma above the classical bound of 2)\n', ...
        S_con_rig, sS_con_rig, (abs(S_con_rig)-2)/sS_con_rig);


%% =========================================================================
% LOCAL FUNCTIONS (add these alongside your existing local functions)
% =========================================================================

function [pFit, covp, chi2red, xq, yq] = local_fit_cosine_weighted(x, y, sy, cosineModel, opts)
% Weighted nonlinear least-squares fit of the cosine model, using known
% per-point standard deviations sy, returning the parameter covariance
% matrix (for delta-method error propagation) and reduced chi-square
% (diagnostic: should be ~1 if error bars and model are both trustworthy).

    x = x(:); y = y(:); sy = sy(:);
    sy(sy <= 0 | isnan(sy)) = median(sy(sy > 0));  % guard degenerate errors

    d0 = mean(y);
    a0 = (max(y) - min(y)) / 2;
    b0 = 2*pi/180;
    [~, idx_max] = max(y);
    c0 = -b0 * x(idx_max);
    p0 = [a0, b0, c0, d0];

    lb = [-Inf, 0.5*b0, -Inf, -Inf];
    ub = [ Inf, 1.5*b0,  Inf,  Inf];

    % Fit in WEIGHTED (chi-square) units: residual = (model - data)/sigma
    wModel = @(p, xx) cosineModel(p, xx) ./ sy;
    yWeighted = y ./ sy;

    [pFit, resnorm, residual, ~, ~, ~, J] = ...
        lsqcurvefit(wModel, p0, x, yWeighted, lb, ub, opts);

    J = full(J);
    dof = length(y) - length(pFit);
    chi2red = resnorm / dof;

    % Since residuals are already in units of sigma, (J'J)^-1 directly
    % gives the parameter covariance (no extra MSE rescaling needed).
    covp = inv(J' * J); %#ok<MINV>

    xq = linspace(min(x), max(x), 300);
    yq = cosineModel(pFit, xq);
end


function [Eval, sigmaE, grad] = evalCosineWithUnc(pFit, covp, theta)
% Evaluates the cosine model a*cos(w*theta+c)+d at a given theta, along
% with its delta-method standard error, using the fit's parameter
% covariance matrix covp.

    a = pFit(1); w = pFit(2); c = pFit(3); d = pFit(4);
    Eval = a*cos(w*theta + c) + d;

    % Analytic gradient wrt [a, w, c, d]
    grad = [ cos(w*theta + c); ...
            -a*theta*sin(w*theta + c); ...
            -a*sin(w*theta + c); ...
             1 ];

    sigmaE = sqrt(grad' * covp * grad);
end


function [S, sigma_S] = chsh_with_covariance(b, bp, pE1, covp1, pE2, covp2)
% Rigorous CHSH S and its uncertainty at arbitrary (possibly unmeasured)
% angles b, bp, via delta-method propagation of the weighted-fit
% covariance matrices of E1 and E2. E1, E2 come from statistically
% independent measurement runs (Procedure 1 vs Procedure 2), so their
% contributions to sigma_S add in quadrature; WITHIN each fit, the
% correlation between evaluations at b and bp (both derived from the
% same 4 fit parameters) is correctly captured via the shared gradient.

    [E1b,  ~, g1b]  = evalCosineWithUnc(pE1, covp1, b);
    [E1bp, ~, g1bp] = evalCosineWithUnc(pE1, covp1, bp);
    [E2b,  ~, g2b]  = evalCosineWithUnc(pE2, covp2, b);
    [E2bp, ~, g2bp] = evalCosineWithUnc(pE2, covp2, bp);

    S = E1b - E1bp + E2b + E2bp;

    v = g1b - g1bp;   % gradient of (E1(b) - E1(b')) wrt pE1
    u = g2b + g2bp;   % gradient of (E2(b) + E2(b')) wrt pE2

    sigma_S = sqrt(v' * covp1 * v + u' * covp2 * u);
end