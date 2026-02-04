function plot_pareto2d_projection(cfg, results, p_star, outbase)
%PLOT_PARETO2D_PROJECTION 2D Pareto projection plot.
%
% x = J_rank (higher better)
% y = J_sep  (higher better)
% color = D_conv (lower better)
% Highlight:
%   - Pareto front points (under 3 objectives)
%   - Selected p_star

results = results(:);

% ---- collect values ----
pvals = [results.p];
pvals = pvals(:);

Jrank = arrayfun(@(s) s.metrics.J_rank, results);
Jsep  = arrayfun(@(s) s.metrics.J_sep,  results);
Dconv = arrayfun(@(s) s.metrics.D_conv, results);
Dconv(isnan(Dconv)) = 0;

% Pareto check: maximize [Jrank, Jsep, -Dconv]
A = [Jrank(:), Jsep(:), -Dconv(:)];
isPareto = pareto_front_maximize(A, 1e-12);

idx_star = find(pvals == p_star, 1);
if isempty(idx_star)
    warning('p_star=%d not found in results. Will not mark star.', p_star);
end

% ---- plot ----
fig = figure('Visible','off');
plot_theme(cfg, fig, true);
hold on;

% base scatter (color by Dconv)
hAll = scatter(Jrank(:), Jsep(:), 90, Dconv(:), 'filled');
hAll.MarkerEdgeColor = [0.20 0.20 0.20];
hAll.LineWidth = 0.8;

colormap(parula);
cb = colorbar;
cb.Label.String = 'D_{conv} (lower better)';

% Pareto ring overlay
scatter(Jrank(isPareto), Jsep(isPareto), 150, 'o', ...
    'MarkerEdgeColor', [0.20 0.45 0.85], 'LineWidth', 2.2);

% p_star marker
if ~isempty(idx_star)
    scatter(Jrank(idx_star), Jsep(idx_star), 220, 'p', ...
        'MarkerFaceColor', [0.85 0.15 0.15], ...
        'MarkerEdgeColor', [0.85 0.15 0.15]);
end

% label only Pareto + p_star
toLabel = isPareto;
if ~isempty(idx_star), toLabel(idx_star) = true; end

for i = 1:numel(pvals)
    if ~toLabel(i), continue; end
    text(Jrank(i), Jsep(i), sprintf('  p=%d', pvals(i)), ...
        'FontSize', 12, 'FontWeight','bold', 'Color', [0 0 0], 'Clipping','on');
end

xlabel('J_{rank}  (higher better)');
ylabel('J_{sep}   (higher better)');
title('2D Pareto Projection (color = D_{conv})', 'FontSize', cfg.plot.title_size);

grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);
box on;

% export
exportgraphics(fig, [outbase '.png'], 'Resolution', cfg.viz.export_dpi);
if cfg.viz.save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType', 'vector');
end
close(fig);

end

% ---------- helper: pareto front for maximize-all ----------
function isPareto = pareto_front_maximize(A, eps)
% A: N x M, maximize all M objectives
N = size(A,1);
isPareto = true(N,1);

for i = 1:N
    if ~isPareto(i), continue; end
    for j = 1:N
        if i == j, continue; end
        ge_all = all(A(j,:) >= A(i,:) - eps);
        gt_one = any(A(j,:) >  A(i,:) + eps);
        if ge_all && gt_one
            isPareto(i) = false;
            break;
        end
    end
end
end
