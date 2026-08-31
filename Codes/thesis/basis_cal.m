% =========================================================================
% TFG - BASIS CALIBRATION / REFERENCE-STATE VERIFICATION: RESULTS PIPELINE
% =========================================================================
% ÍNDEX:
%   SECTION 1: Settings & color palette
%   SECTION 2: File discovery & per-config processing
%   SECTION 3: Diagnostic printout (dominant channel vs. filename)
%   SECTION 4: Summary table (normalized fractions, both error estimates)
%   SECTION 5: Figure - grid of bar charts, one per (pump, basis) config
% =========================================================================

%% =========================================================================
% SECTION 1: Settings & color palette
% =========================================================================
clear; close all; format long g
set(groot, 'defaultTextInterpreter',         'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',       'latex');

tau = 400 * 1e-12;   % [s] coincidence window (from file header)
T_bin = 1;            % [s] each row is a 1-second integration
nBins_expected = 60;  % [PLACEHOLDER: confirm every file has 60 rows]

% Reusing the TE = deep blue / TM = maroon palette from the CAR figures
c_TM = [0.10, 0.20, 0.55];
c_TE = [0.55, 0.10, 0.15];

% [PLACEHOLDER: confirm pump-index convention - assumed 0 = TE, 1 = TM]
pumpIdxLabel = containers.Map({'1','0'}, {'TE', 'TM'});

data_path = '../../Tests/timetagger/wvg2/entr3/';
outFolder = 'results_figures';
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

%% =========================================================================
% SECTION 2: File discovery & per-config processing
% =========================================================================
files = dir(fullfile(data_path, '*.txt'));
fprintf('Found %d files in %s\n', length(files), data_path);

% Group files by (pumpIdx, basisLabel) - e.g. '0_HH', '1_DD', etc.
configs = containers.Map();

for i = 1:length(files)
    fname = files(i).name;
    tok = regexp(fname, '-([A-Z]{2})(\d+)\.txt$', 'tokens');
    if isempty(tok)
        warning('Filename does not match expected pattern, skipping: %s', fname);
        continue
    end
    basisLabel = tok{1}{1};   % e.g. 'HH', 'VV', 'DD', 'AA'
    pumpIdx    = tok{1}{2};   % e.g. '0', '1'
    key = [pumpIdx '_' basisLabel];

    if ~isKey(configs, key)
        configs(key) = {};
    end
    tmp = configs(key);
    tmp{end+1} = fullfile(files(i).folder, fname); %#ok<AGROW>
    configs(key) = tmp;
end

configKeys = keys(configs);
nConfigs = length(configKeys);
fprintf('Discovered %d distinct (pump, basis) configurations:\n', nConfigs);
disp(configKeys');

% --- Process each configuration ---
results = struct([]);

for k = 1:nConfigs
    key = configKeys{k};
    parts = strsplit(key, '_');
    pumpIdx = parts{1};
    basisLabel = parts{2};

    filePaths = configs(key);
    R = processBasisConfig(filePaths, tau);
    R.key = key;
    R.pumpIdx = pumpIdx;
    R.pumpLabel = pumpIdxLabel(pumpIdx);
    R.basisLabel = basisLabel;
    R.nFiles = length(filePaths);

    results = [results, R]; %#ok<AGROW>
end

%% =========================================================================
% SECTION 3: Diagnostic printout (dominant channel vs. filename)
% =========================================================================
fprintf('\n=== DIAGNOSTIC: dominant channel per configuration ===\n');
fprintf('(H/V-setting files SHOULD show one dominant channel near 1;\n');
fprintf(' D/A-setting files SHOULD show all four channels near 1/4.)\n\n');

channelNames = {'HH', 'VV', 'HV', 'VH'};

for k = 1:length(results)
    R = results(k);
    [maxFrac, maxIdx] = max(R.fraction_mean);
    fprintf('%-8s (%d file(s), %d total bins): dominant = %s (%.1f%%) | ', ...
        R.key, R.nFiles, R.nBinsTotal, channelNames{maxIdx}, 100*maxFrac);
    fprintf('all fractions: HH=%.2f VV=%.2f HV=%.2f VH=%.2f\n', R.fraction_mean);

    % Flag mismatch for H/V-type settings (HH or VV filename)
    if any(strcmp(R.basisLabel, {'HH','VV'}))
        expectedChannel = find(strcmp(channelNames, R.basisLabel));
        if maxIdx ~= expectedChannel
            warning(['%s: dominant channel is %s, NOT %s as the filename ' ...
                     'would suggest. This needs investigation before ' ...
                     'writing up Results (see chat discussion).'], ...
                     R.key, channelNames{maxIdx}, R.basisLabel);
        end
    end
end

%% =========================================================================
% SECTION 4: Summary table (normalized fractions, both error estimates)
% =========================================================================
PumpPol = {}; Basis = {}; Channel = {}; MeanTrueRate_cps = [];
Fraction = []; Fraction_errPoisson = []; Fraction_errEmpirical = [];

for k = 1:length(results)
    R = results(k);
    for c = 1:4
        PumpPol{end+1}              = R.pumpLabel; %#ok<AGROW>
        Basis{end+1}                = R.basisLabel; %#ok<AGROW>
        Channel{end+1}              = channelNames{c}; %#ok<AGROW>
        MeanTrueRate_cps(end+1)     = R.trueRate_mean(c); %#ok<AGROW>
        Fraction(end+1)             = R.fraction_mean(c); %#ok<AGROW>
        Fraction_errPoisson(end+1)  = R.fraction_errPoisson(c); %#ok<AGROW>
        Fraction_errEmpirical(end+1)= R.fraction_errEmpirical(c); %#ok<AGROW>
    end
end

T = table(PumpPol', Basis', Channel', MeanTrueRate_cps', Fraction', ...
          Fraction_errPoisson', Fraction_errEmpirical', ...
    'VariableNames', {'PumpPolarization','AnalyzerSetting','Channel', ...
                       'MeanTrueRate_cps','Fraction', ...
                       'Fraction_err_Poisson','Fraction_err_Empirical'});

disp('=== Basis calibration summary table ===');
disp(T);
writetable(T, fullfile(outFolder, 'BasisCalibration_summary.csv'));

%% =========================================================================
% SECTION 5: Figure - grid of bar charts, one per (pump, basis) config
% =========================================================================
% Order panels: TE-H, TE-V, TE-D, TE-A, TM-H, TM-V, TM-D, TM-A (as available)
% basisOrder = {'HH','VV','DD','AA'};
% pumpOrder  = {'0','1'};

basisOrder = {'HH','VV','DD','AA'};
pumpOrder  = {'1','0'};

nRows = length(pumpOrder);
nCols = length(basisOrder);

fig = figure('Name', 'Basis calibration - normalized fractions', ...
             'Color', 'w', 'Position', [60 60 1400 600]);

panelIdx = 1;
for pRow = 1:nRows
    pumpIdx = pumpOrder{pRow};
    for bCol = 1:nCols
        basisLabel = basisOrder{bCol};
        key = [pumpIdx '_' basisLabel];

        ax = subplot(nRows, nCols, panelIdx);
        panelIdx = panelIdx + 1;

        if ~isKey(configs, key)
            title(ax, sprintf('%s - N/A', key), 'FontSize', 10);
            axis(ax, 'off');
            continue
        end

        R = results(strcmp({results.key}, key));
        barColor = c_TE; if strcmp(pumpIdx, '1'), barColor = c_TM; end

        hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        b = bar(ax, 1:4, R.fraction_mean, 'FaceColor', barColor, ...
                'EdgeColor', 'k', 'LineWidth', 0.8);
        errorbar(ax, 1:4, R.fraction_mean, R.fraction_errEmpirical, ...
                 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 6);
        yline(ax, 0.25, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

        set(ax, 'XTick', 1:4, 'XTickLabel', channelNames, 'FontSize', 9);
        ylim(ax, [0 1]);
        title(ax, sprintf('%s pump, %s setting', R.pumpLabel, basisLabel), ...
              'FontSize', 10);
        if bCol == 1, ylabel(ax, 'Fraction', 'FontSize', 10); end
        hold(ax, 'off');
    end
end

%sgtitle('Basis calibration: normalized coincidence fractions per configuration', ...
%        'FontSize', 13, 'Interpreter', 'latex');

saveas(fig, fullfile(outFolder, 'BasisCalibration_Grid.png'));
fprintf('\nAll outputs saved to: %s\n', outFolder);


%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function R = processBasisConfig(filePaths, tau)
% Loads one or more files belonging to the same (pump, basis)
% configuration, computes per-bin true coincidences for all 4 channels,
% and returns mean rates, normalized fractions, and two error estimates
% (Poisson and empirical) for each.

    allTrueHH = []; allTrueVV = []; allTrueHV = []; allTrueVH = [];

    for f = 1:length(filePaths)
        data = readtable(filePaths{f}, 'FileType', 'text', 'Delimiter', ',', ...
            'NumHeaderLines', 11, 'VariableNamingRule', 'preserve');

        S_HA = data{:,3}; S_VA = data{:,4};
        S_HB = data{:,5}; S_VB = data{:,6};

        raw_HH = data{:,8};  raw_VV = data{:,9};
        raw_HV = data{:,10}; raw_VH = data{:,11};

        acc_HH = S_HA .* S_HB .* tau;
        acc_VV = S_VA .* S_VB .* tau;
        acc_HV = S_HA .* S_VB .* tau;
        acc_VH = S_VA .* S_HB .* tau;

        allTrueHH = [allTrueHH; raw_HH - acc_HH]; %#ok<AGROW>
        allTrueVV = [allTrueVV; raw_VV - acc_VV]; %#ok<AGROW>
        allTrueHV = [allTrueHV; raw_HV - acc_HV]; %#ok<AGROW>
        allTrueVH = [allTrueVH; raw_VH - acc_VH]; %#ok<AGROW>
    end

    nBinsTotal = length(allTrueHH);
    T_total = nBinsTotal * 1;   % 1 s per bin

    trueRate_mean = [mean(allTrueHH), mean(allTrueVV), mean(allTrueHV), mean(allTrueVH)];
    trueRate_mean = max(trueRate_mean, 0);   % guard against negative from noise

    % --- Poisson error (supervisor's formula, generalized to T_total) ---
    errPoisson = sqrt(trueRate_mean ./ T_total);

    % --- Empirical error (spread across bins) ---
    errEmpirical = [std(allTrueHH), std(allTrueVV), std(allTrueHV), std(allTrueVH)] / sqrt(nBinsTotal);

    % --- Normalized fractions ---
    total = sum(trueRate_mean);
    fraction_mean = trueRate_mean / total;
    fraction_errPoisson  = errPoisson  / total;
    fraction_errEmpirical= errEmpirical/ total;

    R.trueRate_mean       = trueRate_mean;
    R.fraction_mean       = fraction_mean;
    R.fraction_errPoisson = fraction_errPoisson;
    R.fraction_errEmpirical = fraction_errEmpirical;
    R.nBinsTotal = nBinsTotal;
end