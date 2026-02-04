function main_select_p()
%MAIN_SELECT_P Run Morris best-p search (deterministic) with Successive Halving.
%
% Outputs:
% - p_star in console
% - saves results under results/data/
%
% NOTE:
% - Deterministic model (no noise)
% - Successive Halving for candidate p search
% - Final evaluation on finalists
% - Visualization dataset + plotting (optional via cfg.viz.enable)

% ---- add project paths (src, subfolders) ----
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
addpath(thisDir);
addpath(genpath(fullfile(thisDir, 'src')));

cfg = config();

% Create output dirs
if ~exist(cfg.out.results_dir, 'dir'), mkdir(cfg.out.results_dir); end
if ~exist(cfg.out.data_dir, 'dir'), mkdir(cfg.out.data_dir); end
if isfield(cfg.out, 'fig_dir')
    if ~exist(cfg.out.fig_dir, 'dir'), mkdir(cfg.out.fig_dir); end
end

rng(cfg.seed.master);

% Build deterministic benchmark coefficients
beta = build_beta(cfg.k, cfg.seed.beta, cfg.groups);

% -------------------- Successive Halving search --------------------
[finalists, search_log] = successive_halving_search(cfg, beta);

% -------------------- Final evaluation on finalists ----------------
final_results = struct();
for i = 1:numel(finalists)
    p = finalists(i);
    fprintf('[FINAL] Evaluating p=%d with h=%d, R=%d...\n', p, cfg.final.h, cfg.final.R);
    run_seed_base = cfg.seed.final + 1000*p;
    [metrics, per_run] = compute_metrics_multirun(cfg, beta, p, cfg.final.h, cfg.final.R, run_seed_base);
    final_results(i).p = p; %#ok<AGROW>
    final_results(i).metrics = metrics; %#ok<AGROW>
    final_results(i).per_run = per_run; %#ok<AGROW>
end

% ---- Fill D_conv across FINALISTS (sorted by p) ----
% Reason: compute_metrics_multirun cannot compute D_conv alone (needs adjacent p).
% We compute it here so choose_p_star can use the full constraint rule.
if numel(finalists) >= 2
    [p_sortedF, ~] = sort(finalists(:)', 'ascend');

    med_mu_starF = cell(1, numel(p_sortedF));
    for ii = 1:numel(p_sortedF)
        p = p_sortedF(ii);
        idx = find([final_results.p] == p, 1);
        med_mu_starF{ii} = median(final_results(idx).per_run.mu_star, 1);
    end

    DconvF = nan(1, numel(p_sortedF));
    for ii = 1:numel(p_sortedF)-1
        v = med_mu_starF{ii};
        vnext = med_mu_starF{ii+1};
        rel = abs(v - vnext) ./ (abs(vnext) + cfg.select.eps);
        DconvF(ii) = max(rel);
    end
    DconvF(end) = 0; % last one: treat as converged within finalist set

    for ii = 1:numel(p_sortedF)
        p = p_sortedF(ii);
        idx = find([final_results.p] == p, 1);
        final_results(idx).metrics.D_conv = DconvF(ii);
    end
else
    % Only one finalist: define D_conv = 0
    final_results(1).metrics.D_conv = 0;
end

% Choose p_star using the deterministic multi-criteria rule (on finalists set)
p_star = choose_p_star(cfg, final_results);

% Print summary
fprintf('\n==================== SUMMARY ====================\n');
fprintf('Finalists: %s\n', mat2str(finalists));
for i = 1:numel(final_results)
    p = final_results(i).p;
    m = final_results(i).metrics;
    fprintf('p=%2d | J_rank=%.4f | D_conv=%.4g | J_sep=%.4g\n', p, m.J_rank, m.D_conv, m.J_sep);
end
fprintf('=> p_star = %d\n', p_star);
fprintf('=================================================\n');

% -------------------- Save core data --------------------
save(fullfile(cfg.out.data_dir, 'beta.mat'), 'beta', 'cfg');
save(fullfile(cfg.out.data_dir, 'search_log.mat'), 'search_log');
save(fullfile(cfg.out.data_dir, 'final_results.mat'), 'final_results', 'p_star');

% Robust write for finalists.csv (Windows may lock files opened by Excel/VSCode)
finalists_path = fullfile(cfg.out.data_dir, 'finalists.csv');

% Try to remove existing file (if locked, delete will fail but we handle below)
if exist(finalists_path, 'file')
    try
        delete(finalists_path);
    catch
        % ignore; may be locked
    end
end

% First attempt: write to target path
ok = false;
try
    fid = fopen(finalists_path, 'w');
    if fid ~= -1
        fprintf(fid, 'p\n');
        fprintf(fid, '%d\n', finalists(:));
        fclose(fid);
        ok = true;
    end
catch
    ok = false;
    try, if exist('fid','var') && fid~=-1, fclose(fid); end, catch, end %#ok<CTCH>
end

% Fallback: write to a user-writable temp directory
if ~ok
    fallback_dir = fullfile(tempdir, 'morris_select_p_outputs');
    if ~exist(fallback_dir, 'dir'), mkdir(fallback_dir); end
    fallback_path = fullfile(fallback_dir, 'finalists.csv');

    fid = fopen(fallback_path, 'w');
    if fid == -1
        error('Failed to write finalists.csv to both target path and tempdir(). Please close any program locking the file and check folder permissions.');
    end
    fprintf(fid, 'p\n');
    fprintf(fid, '%d\n', finalists(:));
    fclose(fid);

    warning('Permission denied writing %s. Wrote finalists.csv to fallback path: %s', finalists_path, fallback_path);
else
    fprintf('[SAVE] finalists.csv saved to: %s\n', finalists_path);
end



% Save a compact CSV for final metrics (with h,R)
% ---- Safe write: final_metrics.csv ----
csv_path = fullfile(cfg.out.data_dir, 'final_metrics.csv');
[fid, actual_csv_path] = safe_fopen(csv_path, 'w', 'morris_select_p_outputs');
fprintf(fid, 'p,J_rank,D_conv,J_sep,h,R\n');
for i = 1:numel(final_results)
    p = final_results(i).p;
    m = final_results(i).metrics;
    fprintf(fid, '%d,%.8f,%.8g,%.8g,%d,%d\n', p, m.J_rank, m.D_conv, m.J_sep, m.h, m.R);
end
fclose(fid);
fprintf('[SAVE] final_metrics.csv saved to: %s\n', actual_csv_path);


% -------------------- Visualization dataset + plotting --------------------
if isfield(cfg, 'viz') && cfg.viz.enable
    if ~exist(cfg.out.fig_dir, 'dir'), mkdir(cfg.out.fig_dir); end

    fprintf('\n[VIZ] Building uniform dataset for ALL p with h=%d, R=%d...\n', cfg.viz.h, cfg.viz.R);

    allP = cfg.P(:)';
    viz_results = struct();
    for i = 1:numel(allP)
        p = allP(i);
        run_seed_base = cfg.seed.final + 900000 + 1000*p; % separate seed stream
        [metrics, per_run] = compute_metrics_multirun(cfg, beta, p, cfg.viz.h, cfg.viz.R, run_seed_base);
        viz_results(i).p = p; %#ok<AGROW>
        viz_results(i).metrics = metrics; %#ok<AGROW>
        viz_results(i).per_run = per_run; %#ok<AGROW>
    end

    % Fill D_conv across ALL p (sorted by p)
    [p_sorted, ord] = sort(allP, 'ascend');
    med_mu_star = cell(1, numel(p_sorted));
    for ii = 1:numel(p_sorted)
        med_mu_star{ii} = viz_results(ord(ii)).per_run.med_mu_star;
    end
    Dconv_all = nan(1, numel(p_sorted));
    for ii = 1:numel(p_sorted)-1
        v = med_mu_star{ii};
        vnext = med_mu_star{ii+1};
        rel = abs(v - vnext) ./ (abs(vnext) + cfg.select.eps);
        Dconv_all(ii) = max(rel);
    end
    Dconv_all(end) = 0; % make convergence plot complete

    for ii = 1:numel(p_sorted)
        p = p_sorted(ii);
        idx = find([viz_results.p] == p, 1);
        viz_results(idx).metrics.D_conv = Dconv_all(ii);
    end

    % Save full viz dataset
    allP_used = allP; %#ok<NASGU>
    save(fullfile(cfg.out.data_dir, 'viz_results_allp.mat'), 'viz_results', 'allP_used');

    % Save summary CSV (all p)
    % ---- Safe write: viz_metrics_allp.csv ----
    csv_path2 = fullfile(cfg.out.data_dir, 'viz_metrics_allp.csv');
    [fid2, actual_csv_path2] = safe_fopen(csv_path2, 'w', 'morris_select_p_outputs');
    fprintf(fid2, 'p,J_rank,D_conv,J_sep,h,R\n');
    for i = 1:numel(viz_results)
        p = viz_results(i).p;
        m = viz_results(i).metrics;
        fprintf(fid2, '%d,%.8f,%.8g,%.8g,%d,%d\n', p, m.J_rank, m.D_conv, m.J_sep, m.h, m.R);
    end
    fclose(fid2);
    fprintf('[SAVE] viz_metrics_allp.csv saved to: %s\n', actual_csv_path2);


    % ---- Plotting ----
    % Functions should exist under src/plotting/
    plot_indices_vs_j(cfg, viz_results, fullfile(cfg.out.fig_dir, 'fig_indices_vs_j'));
    plot_mustar_sigma_multiples(cfg, viz_results, fullfile(cfg.out.fig_dir, 'fig_mustar_sigma'));
    plot_rank_stability(cfg, viz_results, fullfile(cfg.out.fig_dir, 'fig_rank_stability_vs_p'));
    plot_convergence(cfg, viz_results, fullfile(cfg.out.fig_dir, 'fig_convergence_vs_p'));
    plot_pareto(cfg, viz_results, fullfile(cfg.out.fig_dir, 'fig_pareto_cost_perf'));

    selected_factors = [2, 9, 15];
    plot_EE_boxplot(cfg, beta, selected_factors, fullfile(cfg.out.fig_dir, 'fig_EE_boxplot_selected_factors'));
    plot_pareto3d(cfg, viz_results, p_star, fullfile(cfg.out.fig_dir, 'fig_pareto_3d_objectives'));
    plot_pareto2d_projection(cfg, viz_results, p_star, fullfile(cfg.out.fig_dir, 'fig_pareto_2d_projection'));
    plot_muabs_grid_24(cfg, viz_results, fullfile(cfg.out.fig_dir, 'fig_muabs_grid_24'));



    fprintf('[VIZ] Figures saved to: %s\n', cfg.out.fig_dir);
end
% -------------------------------------------------------------------------

end
