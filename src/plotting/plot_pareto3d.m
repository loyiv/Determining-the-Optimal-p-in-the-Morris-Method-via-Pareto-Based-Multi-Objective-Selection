function plot_pareto3d(cfg, results, p_star, outbase)
%PLOT_PARETO3D 3D Pareto front visualization for selecting p.
%
% Objectives:
%   maximize J_rank
%   maximize J_sep
%   minimize D_conv
%
% For dominance checks, we transform to "all maximize" by using Z = -D_conv.

% Ensure column struct array
results = results(:);

% ---- collect objective values ----
pvals = [results.p];
pvals = pvals(:);

Jrank = arrayfun(@(s) s.metrics.J_rank, results);
Jsep  = arrayfun(@(s) s.metrics.J_sep,  results);
Dconv = arrayfun(@(s) s.metrics.D_conv, results);

% sanitize NaN
Dconv(isnan(Dconv)) = 0;

% transform to maximize-all
X = Jrank(:);
Y = Jsep(:);
Z = -Dconv(:);

% Pareto front (maximize-all)
isPareto = pareto_front_maximize([X Y Z], 1e-12);

% locate p_star
idx_star = find(pvals == p_star, 1);
if isempty(idx_star)
    warning('p_star=%d not found in results. Will not mark star.', p_star);
end

% ---- plot ----
fig = figure('Visible','off');
plot_theme(cfg, fig, true);

hold on;

% all points (gray)
scatter3(X, Y, Z, 85, [0.65 0.65 0.65], 'filled');

% pareto points (blue)
scatter3(X(isPareto), Y(isPareto), Z(isPareto), 125, [0.20 0.45 0.85], 'filled');

% p_star marker
if ~isempty(idx_star)
    scatter3(X(idx_star), Y(idx_star), Z(idx_star), 220, [0.85 0.15 0.15], 'p', 'filled');
end

% label only Pareto points + p_star
toLabel = isPareto;
if ~isempty(idx_star), toLabel(idx_star) = true; end

for i = 1:numel(pvals)
    if ~toLabel(i), continue; end
    txt = sprintf('  p=%d', pvals(i));
    text(X(i), Y(i), Z(i), txt, 'FontSize', 12, 'Color', [0 0 0], ...
        'FontWeight','bold', 'Clipping', 'on');
end

xlabel('J_{rank}  (higher better)');
ylabel('J_{sep}   (higher better)');
zlabel('-D_{conv} (higher better)');
title('3D Pareto Front in Objective Space', 'FontSize', cfg.plot.title_size);

grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);
view(45, 22);
box on;

% legend (proxy)
h1 = plot(nan,nan,'o','MarkerSize',9,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none');
h2 = plot(nan,nan,'o','MarkerSize',9,'MarkerFaceColor',[0.20 0.45 0.85],'MarkerEdgeColor','none');
h3 = plot(nan,nan,'p','MarkerSize',12,'MarkerFaceColor',[0.85 0.15 0.15],'MarkerEdgeColor','none');
lg = legend([h1 h2 h3], {'All candidates','Pareto front','Selected p^*'}, 'Location','best');
lg.Box = 'off';

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
