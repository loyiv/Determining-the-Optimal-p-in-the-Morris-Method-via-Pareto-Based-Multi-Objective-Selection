function export_EE_long_csv()
%EXPORT_EE_LONG_CSV Export raw Morris elementary effects (EE) into a long-form CSV.
%
% Output CSV columns:
%   p, factor, traj, ee
%
% This script REUSES your existing implementation:
%   [EE, stats] = morris_run(cfg, beta, p, h, seed)
%
% Requirements:
%   - config.m provides cfg with fields:
%       cfg.k, cfg.P, cfg.viz.h (or cfg.ee_export.h_override),
%       cfg.out.results_dir, cfg.out.data_dir, cfg.seed.final
%   - beta.mat exists under cfg.out.data_dir and contains variable beta
%   - morris_run.m is on MATLAB path (typically under src/)
%
% Output directory:
%   <cfg.out.results_dir>/csv_export/fig_EE_boxplot_selected_factors/ee_long.csv
% If permission denied, falls back to:
%   tempdir()/morris_select_p_outputs/csv_export/fig_EE_boxplot_selected_factors/ee_long.csv

% ---- project paths ----
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
addpath(thisDir);
addpath(genpath(fullfile(thisDir, 'src')));

cfg = config();

% ---- load beta ----
betaMat = fullfile(cfg.out.data_dir, 'beta.mat');
if ~exist(betaMat, 'file')
    error('Cannot find beta.mat at: %s. Please run main_select_p first.', betaMat);
end
S = load(betaMat);
if ~isfield(S, 'beta')
    error('beta.mat does not contain variable "beta".');
end
beta = S.beta;

% ---- choose p set ----
allP = cfg.P(:)';  % default: all candidates
if isfield(cfg, 'ee_export') && isfield(cfg.ee_export, 'P_subset') && ~isempty(cfg.ee_export.P_subset)
    allP = cfg.ee_export.P_subset(:)';
end

% ---- selected factors ----
if isfield(cfg, 'ee_export') && isfield(cfg.ee_export, 'selected_factors') && ~isempty(cfg.ee_export.selected_factors)
    selected_factors = cfg.ee_export.selected_factors(:)';
else
    % default consistent with your previous choice
    selected_factors = [2, 9, 15];
end

% ---- trajectories per p for EE export ----
h = cfg.viz.h;
if isfield(cfg, 'ee_export') && isfield(cfg.ee_export, 'h_override') && ~isempty(cfg.ee_export.h_override)
    h = cfg.ee_export.h_override;
end

% ---- seed base ----
baseSeed = cfg.seed.final + 700000;

% ---- output dir ----
exportDirPref = fullfile(cfg.out.results_dir, 'csv_export', 'fig_EE_boxplot_selected_factors');
exportDir = ensure_writable_dir(exportDirPref);
csvPathPref = fullfile(exportDir, 'ee_long.csv');

fprintf('[EE-EXPORT] Exporting raw EE samples...\n');
fprintf('  P = %s\n', mat2str(allP));
fprintf('  selected_factors = %s\n', mat2str(selected_factors));
fprintf('  h (trajectories per p) = %d\n', h);
fprintf('  preferred output = %s\n', csvPathPref);

% ---- open file safely (permission fallback) ----
[fid, actualPath] = safe_fopen(csvPathPref, 'w', 'morris_select_p_outputs');
fprintf(fid, 'p,factor,traj,ee\n');

% ---- main loop ----
for ip = 1:numel(allP)
    p = allP(ip);
    seed = baseSeed + 1000*p;

    % Reuse your existing Morris implementation:
    % Expect EE to be either k x h or h x k
    [EE_all, ~] = morris_run(cfg, beta, p, h, seed);

    % adapt shape -> k x h
    if size(EE_all, 1) == h && size(EE_all, 2) == cfg.k
        EE_all = EE_all.'; % h x k -> k x h
    end

    if size(EE_all, 1) ~= cfg.k
        fclose(fid);
        error('EE_all first dimension must be k=%d after adaptation, got %d (p=%d).', cfg.k, size(EE_all,1), p);
    end
    if size(EE_all, 2) ~= h
        fclose(fid);
        error('EE_all second dimension must be h=%d after adaptation, got %d (p=%d).', h, size(EE_all,2), p);
    end

    % write only selected factors
    for jf = 1:numel(selected_factors)
        j = selected_factors(jf);
        if j < 1 || j > cfg.k
            fclose(fid);
            error('Selected factor %d out of range [1,%d].', j, cfg.k);
        end
        ee = EE_all(j, :);

        for r = 1:h
            fprintf(fid, '%d,%d,%d,%.12g\n', p, j, r, ee(r));
        end
    end

    fprintf('[EE-EXPORT] Done p=%d (%d/%d)\n', p, ip, numel(allP));
end

fclose(fid);
fprintf('[EE-EXPORT] Wrote ee_long.csv to: %s\n', actualPath);

end


% ======================= helper: ensure writable dir =======================
function outDir = ensure_writable_dir(preferredDir)
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
    outDir = fullfile(tempdir, 'morris_select_p_outputs', 'csv_export', 'fig_EE_boxplot_selected_factors');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
end
end


% ======================= helper: safe fopen with fallback =======================
function [fid, actual_path] = safe_fopen(preferred_path, mode, fallback_subdir)
if nargin < 3 || isempty(fallback_subdir)
    fallback_subdir = 'morris_select_p_outputs';
end

pref_dir = fileparts(preferred_path);
if ~isempty(pref_dir) && ~exist(pref_dir, 'dir')
    mkdir(pref_dir);
end

fid = fopen(preferred_path, mode);
if fid ~= -1
    actual_path = preferred_path;
    return;
end

fallback_dir = fullfile(tempdir, fallback_subdir, 'csv_export', 'fig_EE_boxplot_selected_factors');
if ~exist(fallback_dir, 'dir'), mkdir(fallback_dir); end

[~, name, ext] = fileparts(preferred_path);
fallback_path = fullfile(fallback_dir, [name ext]);

fid = fopen(fallback_path, mode);
if fid == -1
    error('SAFE_FOPEN:Failed', ...
        'Cannot open file for writing:\nPreferred: %s\nFallback:  %s\n', ...
        preferred_path, fallback_path);
end

actual_path = fallback_path;
warning('Permission denied writing %s. Using fallback: %s', preferred_path, fallback_path);
end
