function plot_EE_boxplot(cfg, beta, factors, outbase)
%PLOT_EE_BOXPLOT Boxplot of EE distributions for selected factors across p.

P = cfg.P(:)';
[~, ord] = sort(P, 'ascend');
P = P(ord);

h = cfg.viz.h;
Delta = zeros(size(P));
for i = 1:numel(P), Delta(i) = 2/(P(i)-1); end %#ok<NASGU>

% For each p, run once with h trajectories to collect EE
EE_all = cell(numel(P),1);
for i = 1:numel(P)
    p = P(i);
    run_seed = cfg.seed.final + 777000 + 1000*p;
    [EE, ~] = morris_run(cfg, beta, p, h, run_seed);
    EE_all{i} = EE;
end

fig = figure('Visible','off');
plot_theme(cfg, fig, false);

t = tiledlayout(numel(factors), 1, 'TileSpacing','compact','Padding','compact');

for fidx = 1:numel(factors)
    j = factors(fidx);

    % stack data for boxplot: groups = p
    data = [];
    group = [];
    for i = 1:numel(P)
        ee = EE_all{i}(:, j);
        data = [data; ee]; %#ok<AGROW>
        group = [group; repmat(P(i), numel(ee), 1)]; %#ok<AGROW>
    end

    nexttile; hold on;
    boxplot(data, group, 'Symbol','.', 'Widths', 0.7);
    title(sprintf('EE distribution for factor j=%d', j), 'FontSize', cfg.plot.title_size);
    xlabel('p'); ylabel('EE');
    grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);
end

exportgraphics(fig, [outbase '.png'], 'Resolution', cfg.viz.export_dpi);
if cfg.viz.save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType','vector');
end
close(fig);

end
