function plot_indices_vs_j(cfg, results, outbase)
%PLOT_INDICES_VS_J 3-panel plot: |mu|, mu_star, sigma vs factor index.
% results: struct array with fields p, per_run.(mu, mu_star, sigma)

addpath(fileparts(mfilename('fullpath')));

P = [results.p];
[~, ord] = sort(P, 'ascend');
results = results(ord);

k = cfg.k;
j = 1:k;

% prepare curves using median across repeats
mu_abs = zeros(numel(results), k);
mu_star = zeros(numel(results), k);
sigma = zeros(numel(results), k);

for i = 1:numel(results)
    mu_abs(i,:) = median(abs(results(i).per_run.mu), 1);
    mu_star(i,:) = median(results(i).per_run.mu_star, 1);
    sigma(i,:) = median(results(i).per_run.sigma, 1);
end

fig = figure('Visible','off');
plot_theme(cfg, fig, false);

t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Colormap for p curves
cm = lines(numel(results));

nexttile; hold on;
for i = 1:numel(results)
    plot(j, mu_abs(i,:), '-', 'Color', cm(i,:));
end
title('| \mu_j(p) | vs factor index', 'FontSize', cfg.plot.title_size);
xlabel('Factor j'); ylabel('|mu|');
grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);

nexttile; hold on;
for i = 1:numel(results)
    plot(j, mu_star(i,:), '-', 'Color', cm(i,:));
end
title('\mu_j^*(p) vs factor index', 'FontSize', cfg.plot.title_size);
xlabel('Factor j'); ylabel('mu*');
grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);

nexttile; hold on;
for i = 1:numel(results)
    plot(j, sigma(i,:), '-', 'Color', cm(i,:));
end
title('\sigma_j(p) vs factor index', 'FontSize', cfg.plot.title_size);
xlabel('Factor j'); ylabel('sigma');
grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);

lg = legend(arrayfun(@(x) sprintf('p=%d', x.p), results, 'UniformOutput', false), ...
    'Location','eastoutside');
lg.Box = 'off';

% export
exportgraphics(fig, [outbase '.png'], 'Resolution', cfg.viz.export_dpi);
if cfg.viz.save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType', 'vector');
end
close(fig);
end
