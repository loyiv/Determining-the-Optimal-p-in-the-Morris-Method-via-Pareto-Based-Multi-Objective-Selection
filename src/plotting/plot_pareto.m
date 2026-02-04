function plot_pareto(cfg, results, outbase)
%PLOT_PARETO Plot cost vs stability (or cost vs composite) for all p.

P = [results.p];
[~, ord] = sort(P, 'ascend');
results = results(ord);

pvals = [results.p];
Jrank = arrayfun(@(s) s.metrics.J_rank, results);
h = arrayfun(@(s) s.metrics.h, results);
R = arrayfun(@(s) s.metrics.R, results);

% deterministic no-noise: cost proportional to R*h*(k+1)
cost = R .* h .* (cfg.k + 1);

fig = figure('Visible','off');
plot_theme(cfg, fig, true);
hold on;

scatter(cost, Jrank, 80, [0.25 0.25 0.25], 'filled');
for i = 1:numel(pvals)
    text(cost(i), Jrank(i), sprintf('  p=%d', pvals(i)), 'FontSize', 11);
end

xlabel('Cost ~ R \cdot h \cdot (k+1)');
ylabel('J_{rank}');
title('Pareto view: stability vs cost', 'FontSize', cfg.plot.title_size);
grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);

exportgraphics(fig, [outbase '.png'], 'Resolution', cfg.viz.export_dpi);
if cfg.viz.save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType','vector');
end
close(fig);
end
