function export_viz_to_csv()
%EXPORT_VIZ_TO_CSV Export all visualization datasets to CSV files (per-figure).
%
% This script reads:
%   results/data/viz_results_allp.mat  (contains viz_results)
%   results/data/beta.mat (optional)   (contains beta, cfg)
%
% Outputs:
%   results/csv_export/...
% If results/csv_export is not writable, falls back to tempdir().
%
% Each "figure" gets its own CSV(s), so you can re-plot in Python.

% -------- paths --------
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
addpath(thisDir);
addpath(genpath(fullfile(thisDir, 'src')));

cfg = config();

% -------- robust locate viz_results_allp.mat --------
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);

% Candidate locations (try in order)
cands = {
    fullfile(cfg.out.data_dir, 'viz_results_allp.mat')                                      % from config()
    fullfile(thisDir, 'results', 'data', 'viz_results_allp.mat')                           % project-local results
    fullfile(fileparts(thisDir), 'results', 'data', 'viz_results_allp.mat')                % one level up
    fullfile(thisDir, 'final', 'results', 'data', 'viz_results_allp.mat')                  % common "final/" nesting
    fullfile(fileparts(thisDir), 'final', 'results', 'data', 'viz_results_allp.mat')       % one level up + final/
    fullfile(tempdir, 'morris_select_p_outputs', 'viz_results_allp.mat')                   % temp fallback (if you ever saved there)
    fullfile(tempdir, 'morris_select_p_outputs', 'results', 'data', 'viz_results_allp.mat')% another temp layout
};

vizMat = '';
for i = 1:numel(cands)
    if exist(cands{i}, 'file')
        vizMat = cands{i};
        break;
    end
end

if isempty(vizMat)
    fprintf(2, 'Cannot find viz_results_allp.mat. Tried:\n');
    for i = 1:numel(cands)
        fprintf(2, '  - %s\n', cands{i});
    end
    error(['Cannot find viz_results_allp.mat. ', ...
           'Please ensure main_select_p ran with cfg.viz.enable=true and saved viz_results_allp.mat.']);
end

fprintf('[EXPORT] Found viz_results_allp.mat at: %s\n', vizMat);

S = load(vizMat);
if ~isfield(S, 'viz_results')
    error('The MAT file does not contain variable "viz_results": %s', vizMat);
end
viz_results = S.viz_results(:);

% Try load beta.mat near vizMat (same folder) for cfg/beta (optional)
dataDirFound = fileparts(vizMat);
betaMat = fullfile(dataDirFound, 'beta.mat');
beta = [];
if exist(betaMat, 'file')
    S2 = load(betaMat);
    if isfield(S2, 'beta'), beta = S2.beta; end
    if isfield(S2, 'cfg'),  cfg  = S2.cfg;  end
end


S = load(vizMat);  % expects viz_results
if ~isfield(S, 'viz_results')
    error('viz_results_allp.mat does not contain variable "viz_results".');
end
viz_results = S.viz_results(:);

% load beta & cfg if exists (for some exports)
beta = [];
if exist(betaMat, 'file')
    S2 = load(betaMat);
    if isfield(S2, 'beta'), beta = S2.beta; end
    if isfield(S2, 'cfg'), cfg = S2.cfg; end
end

% -------- choose export root --------
exportRootPreferred = fullfile(cfg.out.results_dir, 'csv_export');
exportRoot = ensure_writable_dir(exportRootPreferred);

fprintf('[EXPORT] Writing CSV exports to: %s\n', exportRoot);

% sort by p
p_all = [viz_results.p]';
[p_sorted, ord] = sort(p_all, 'ascend');
viz_sorted = viz_results(ord);

% ============================================================
% 0) Global metrics table (already have, but export again)
% ============================================================
dir0 = ensure_dir(fullfile(exportRoot, 'metrics'));
T_metrics = table();
T_metrics.p      = p_sorted;
T_metrics.J_rank = arrayfun(@(s) s.metrics.J_rank, viz_sorted);
T_metrics.D_conv = arrayfun(@(s) s.metrics.D_conv, viz_sorted);
T_metrics.J_sep  = arrayfun(@(s) s.metrics.J_sep,  viz_sorted);
T_metrics.h      = arrayfun(@(s) s.metrics.h,       viz_sorted);
T_metrics.R      = arrayfun(@(s) s.metrics.R,       viz_sorted);

writetable(T_metrics, fullfile(dir0, 'metrics_allp.csv'));

% ============================================================
% 1) Figure: indices vs factor index (3 panels)
% Export long-form and wide-form:
%   - indices_long.csv: columns [p, j, mu, mu_abs, mu_star, sigma]
%   - indices_wide_mu_star.csv etc: each p as a column
% ============================================================
% ============================================================
% 1) Figure: indices vs factor index (3 panels)
% Export long-form and wide-form:
%   - indices_long.csv: columns [p, j, mu, mu_abs, mu_star, sigma]
%     (mu/mu_abs may be NaN if med_mu is not stored in viz_results)
% ============================================================
dir1 = ensure_dir(fullfile(exportRoot, 'fig_indices_vs_j'));

% detect k from med_mu_star length
k = numel(viz_sorted(1).per_run.med_mu_star);

rows = numel(p_sorted) * k;
p_col = zeros(rows,1);
j_col = zeros(rows,1);

mu_col     = nan(rows,1);   % may be unavailable
muabs_col  = nan(rows,1);   % may be unavailable
mustar_col = nan(rows,1);
sig_col    = nan(rows,1);

t = 1;

% check available fields once
prTest = viz_sorted(1).per_run;
has_med_mu    = isfield(prTest, 'med_mu');
has_med_muabs = isfield(prTest, 'med_mu_abs');  % optional alternative

if ~has_med_mu && ~has_med_muabs
    warning(['per_run does not contain med_mu (or med_mu_abs). ', ...
             'Will export mu/mu_abs as NaN. ', ...
             'If you need the |mu| panel in Python, patch compute_metrics_multirun to store med_mu.']);
end

for i = 1:numel(p_sorted)
    pr = viz_sorted(i).per_run;

    mustar = pr.med_mu_star(:);
    sig    = pr.med_sigma(:);

    % optional: mu and/or mu_abs
    mu = nan(k,1);
    muabs = nan(k,1);

    if isfield(pr, 'med_mu')
        mu = pr.med_mu(:);
        muabs = abs(mu);
    elseif isfield(pr, 'med_mu_abs')
        muabs = pr.med_mu_abs(:);
        % keep mu as NaN (unknown sign)
    end

    for j = 1:k
        p_col(t) = p_sorted(i);
        j_col(t) = j;
        mu_col(t) = mu(j);
        muabs_col(t) = muabs(j);
        mustar_col(t) = mustar(j);
        sig_col(t) = sig(j);
        t = t + 1;
    end
end

T_idx = table(p_col, j_col, mu_col, muabs_col, mustar_col, sig_col, ...
    'VariableNames', {'p','j','mu','mu_abs','mu_star','sigma'});
writetable(T_idx, fullfile(dir1, 'indices_long.csv'));

% wide-form: one file per metric (rows=j, cols=p_4, p_5, ...)
W_mustar = make_wide(p_sorted, k, T_idx, 'mu_star');
W_sigma  = make_wide(p_sorted, k, T_idx, 'sigma');
writetable(W_mustar, fullfile(dir1, 'indices_wide_mu_star.csv'));
writetable(W_sigma,  fullfile(dir1, 'indices_wide_sigma.csv'));

% only write mu_abs wide-form if it is not all NaN
if any(~isnan(T_idx.mu_abs))
    W_muabs = make_wide(p_sorted, k, T_idx, 'mu_abs');
    writetable(W_muabs, fullfile(dir1, 'indices_wide_mu_abs.csv'));
end


% ============================================================
% 2) Figure: mu_star-sigma scatter per p
% Export one CSV per p: columns [j, mu_star, sigma]
% ============================================================
dir2 = ensure_dir(fullfile(exportRoot, 'fig_mustar_sigma_per_p'));
for i = 1:numel(p_sorted)
    pr = viz_sorted(i).per_run;
    mustar = pr.med_mu_star(:);
    sig    = pr.med_sigma(:);
    j = (1:k)';

    T = table(j, mustar, sig, 'VariableNames', {'j','mu_star','sigma'});
    writetable(T, fullfile(dir2, sprintf('mustar_sigma_p_%d.csv', p_sorted(i))));
end

% Also export a combined long-form:
% columns [p, j, mu_star, sigma]
p_rep = repelem(p_sorted, k);
j_rep = repmat((1:k)', numel(p_sorted), 1);
mustar_all = T_idx.mu_star;
sigma_all  = T_idx.sigma;
T_ms = table(p_rep, j_rep, mustar_all, sigma_all, 'VariableNames', {'p','j','mu_star','sigma'});
writetable(T_ms, fullfile(dir2, 'mustar_sigma_long.csv'));

% ============================================================
% 3) Figure: rank stability vs p
% If per_run stores tau distribution, export it.
% If not available, export J_rank only (already in metrics).
% Here we export:
%   - rank_stability_summary.csv (p, J_rank)
%   - rank_stability_runs.csv (p, run, tau) if available
% ============================================================
dir3 = ensure_dir(fullfile(exportRoot, 'fig_rank_stability'));

T_rank = table(p_sorted, T_metrics.J_rank, 'VariableNames', {'p','J_rank'});
writetable(T_rank, fullfile(dir3, 'rank_stability_summary.csv'));

% If your compute_metrics_multirun stored tau_runs (vector length R)
has_tau = isfield(viz_sorted(1).per_run, 'tau_runs');
if has_tau
    pr0 = viz_sorted(1).per_run.tau_runs;
    R = numel(pr0);
    p_col = [];
    r_col = [];
    tau_col = [];
    for i = 1:numel(p_sorted)
        tau = viz_sorted(i).per_run.tau_runs(:);
        p_col = [p_col; repmat(p_sorted(i), R, 1)]; %#ok<AGROW>
        r_col = [r_col; (1:R)']; %#ok<AGROW>
        tau_col = [tau_col; tau]; %#ok<AGROW>
    end
    T_tau = table(p_col, r_col, tau_col, 'VariableNames', {'p','run','tau'});
    writetable(T_tau, fullfile(dir3, 'rank_stability_runs.csv'));
end

% ============================================================
% 4) Figure: convergence vs p
% Export convergence curve (p, D_conv, tau_conv)
% ============================================================
dir4 = ensure_dir(fullfile(exportRoot, 'fig_convergence'));

tau = cfg.select.tau_conv;
T_conv = table(p_sorted, T_metrics.D_conv, repmat(tau, numel(p_sorted), 1), ...
    'VariableNames', {'p','D_conv','tau_conv'});
writetable(T_conv, fullfile(dir4, 'convergence_vs_p.csv'));

% ============================================================
% 5) Figure: Pareto (objective space)
% Export:
%   - pareto_objectives.csv (p, J_rank, J_sep, D_conv, minusDconv)
% ============================================================
dir5 = ensure_dir(fullfile(exportRoot, 'fig_pareto'));

minusD = -T_metrics.D_conv;
T_par = table(p_sorted, T_metrics.J_rank, T_metrics.J_sep, T_metrics.D_conv, minusD, ...
    'VariableNames', {'p','J_rank','J_sep','D_conv','minusD_conv'});
writetable(T_par, fullfile(dir5, 'pareto_objectives.csv'));

fprintf('[EXPORT] Done.\n');
fprintf('Tip: Use Python to read CSVs under: %s\n', exportRoot);

end

% ======================== helpers ========================

function outDir = ensure_writable_dir(preferredDir)
% Create preferredDir; if cannot write, fall back to tempdir
outDir = preferredDir;
try
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    testfile = fullfile(outDir, sprintf('.write_test_%d.tmp', randi(1e9)));
    fid = fopen(testfile, 'w');
    if fid == -1, error('cannot open'); end
    fprintf(fid, 'ok');
    fclose(fid);
    delete(testfile);
catch
    outDir = fullfile(tempdir, 'morris_select_p_outputs', 'csv_export');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
end
end

function d = ensure_dir(d)
if ~exist(d, 'dir'), mkdir(d); end
end

function W = make_wide(p_sorted, k, T_idx, fieldname)
% Build wide table: rows=j, columns=p_4, p_5, ...
j = (1:k)';
W = table(j, 'VariableNames', {'j'});
for i = 1:numel(p_sorted)
    p = p_sorted(i);
    sel = (T_idx.p == p);
    v = T_idx.(fieldname)(sel);
    if numel(v) ~= k
        error('make_wide: expected %d rows for p=%d, got %d', k, p, numel(v));
    end
    colname = matlab.lang.makeValidName(sprintf('p_%d', p));
    W.(colname) = v;
end
end
