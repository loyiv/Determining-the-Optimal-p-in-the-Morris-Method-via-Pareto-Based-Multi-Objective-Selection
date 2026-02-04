function plot_muabs_grid_24(cfg, results, outbase)
%PLOT_MUABS_GRID_24  6x4 grid: |mu_j(p)| vs factor index (24 subplots).
% Each subplot corresponds to one p and has ONE curve.
%
% Also exports ONLY the data used for this figure to:
%   <cfg.out.results_dir>/csv_export/fig_indices_muabs_grid/
%   - muabs_grid_long.csv   (p, j, mu_abs)  <-- recommended for Python
%   - muabs_grid_wide.csv   (j, c01..c24)   + muabs_grid_pmap.csv (c01->p)
%
% NOTE: We DO NOT use p as column names (to avoid array2table VariableNames issues).

% ---------------- sort by p ----------------
P = [results.p];
[~, ord] = sort(P, 'ascend');
results = results(ord);
P = P(ord);

k = cfg.k;
j = (1:k)';

% ---------------- compute |mu| curve per p ----------------
mu_abs_all = zeros(numel(results), k);
for i = 1:numel(results)
    mu_abs_all(i,:) = median(abs(results(i).per_run.mu), 1); % per_run.mu is R x k
end

% ---------------- choose 24 p values (uniform) ----------------
N = numel(results);
target = min(24, N);

idx = unique(round(linspace(1, N, target)));
% fill if not enough unique
if numel(idx) < target
    cand = setdiff(1:N, idx, 'stable');
    idx = [idx, cand(1:(target-numel(idx)))];
end
idx = idx(1:target);

% ensure strictly unique & stable
idx = unique(idx, 'stable');
target = numel(idx);

P_sel = P(idx);
mu_abs_sel = mu_abs_all(idx, :);   % target x k

% ---------------- export ONLY used data ----------------
csvDir = fullfile(cfg.out.results_dir, 'csv_export', 'fig_indices_muabs_grid');
if ~exist(csvDir, 'dir'), mkdir(csvDir); end

% Long format (recommended)
[PP, JJ] = ndgrid(P_sel(:), (1:k)');
vals = mu_abs_sel.';   % k x target
Tlong = table(PP(:), JJ(:), vals(:), 'VariableNames', {'p','j','mu_abs'});
writetable(Tlong, fullfile(csvDir, 'muabs_grid_long.csv'));

% Wide format with safe column names c01..cXX
colNames = arrayfun(@(t) sprintf('c%02d', t), 1:target, 'UniformOutput', false);
Twide = array2table(mu_abs_sel.', 'VariableNames', colNames); % (k x target)
Twide = addvars(Twide, j, 'Before', 1, 'NewVariableNames', 'j');
writetable(Twide, fullfile(csvDir, 'muabs_grid_wide.csv'));

% Column -> p mapping
Tpmap = table(colNames(:), P_sel(:), 'VariableNames', {'col','p'});
writetable(Tpmap, fullfile(csvDir, 'muabs_grid_pmap.csv'));

fprintf('[CSV] Saved mu_abs grid data to: %s\n', csvDir);

% ---------------- plot ----------------
title_size = get_field(cfg, {'plot','title_size'}, 12);
axis_size  = get_field(cfg, {'plot','axis_size'}, 10);
dpi        = get_field(cfg, {'viz','export_dpi'}, 300);
save_pdf   = get_field(cfg, {'viz','save_pdf'}, true);

fig = figure('Visible','off');
set(fig, 'Color', 'w');

tlo = tiledlayout(fig, 6, 4, 'TileSpacing','compact', 'Padding','compact');

for ii = 1:target
    ax = nexttile(tlo);
    hold(ax, 'on');

    y = mu_abs_sel(ii,:); % 1 x k
    plot(ax, 1:k, y, '-', 'LineWidth', 2.0, 'Color', [0 0 0]);

    title(ax, sprintf('|\\mu_j(p)| vs factor index  (p=%d)', P_sel(ii)), 'FontSize', title_size);
    xlabel(ax, 'Factor j', 'FontSize', axis_size);
    ylabel(ax, '|mu|', 'FontSize', axis_size);

    grid(ax, 'off');
    ax.LineWidth = 1.2;
    box(ax, 'on');
end

% export
outDir = fileparts(outbase);
if ~exist(outDir, 'dir'), mkdir(outDir); end
exportgraphics(fig, [outbase '.png'], 'Resolution', dpi);
if save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType', 'vector');
end
close(fig);

fprintf('[FIG] Saved: %s(.png/.pdf)\n', outbase);

end

% -------- helper: safe nested field getter --------
function v = get_field(s, path, default)
v = default;
try
    cur = s;
    for i = 1:numel(path)
        key = path{i};
        if ~isstruct(cur) || ~isfield(cur, key), return; end
        cur = cur.(key);
    end
    if ~isempty(cur), v = cur; end
catch
    v = default;
end
end
