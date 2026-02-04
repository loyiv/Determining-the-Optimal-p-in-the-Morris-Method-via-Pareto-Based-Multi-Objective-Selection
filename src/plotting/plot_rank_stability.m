function plot_rank_stability(cfg, results, outbase)
%PLOT_RANK_STABILITY Plot J_rank vs p with IQR error band estimated from pairwise taus.

P = [results.p];
[~, ord] = sort(P, 'ascend');
results = results(ord);
pvals = [results.p];

J = zeros(size(pvals));
iqr_lo = zeros(size(pvals));
iqr_hi = zeros(size(pvals));

for i = 1:numel(results)
    rank_mat = results(i).per_run.rank; % R x k
    R = size(rank_mat,1);

    taus = [];
    idx = 1;
    for a = 1:R-1
        for b = a+1:R
            taus(idx,1) = kendall_tau_b(rank_mat(a,:), rank_mat(b,:)); %#ok<AGROW>
            idx = idx + 1;
        end
    end

    J(i) = median(taus);
    q1 = prctile(taus, 25);
    q3 = prctile(taus, 75);
    iqr_lo(i) = J(i) - q1;
    iqr_hi(i) = q3 - J(i);
end

fig = figure('Visible','off');
plot_theme(cfg, fig, true);
hold on;

% error band (IQR)
x = pvals;
y = J;
y1 = y - iqr_lo;
y2 = y + iqr_hi;
fill([x, fliplr(x)], [y1, fliplr(y2)], [0.85 0.88 0.95], 'EdgeColor','none', 'FaceAlpha',0.6);
plot(x, y, '-o', 'Color', [0.15 0.25 0.55], 'MarkerFaceColor',[0.15 0.25 0.55]);

yline(cfg.select.tau_rank, '--', sprintf('\\tau_{rank}=%.2f', cfg.select.tau_rank), 'LineWidth',1.5);

xlabel('p'); ylabel('J_{rank} (median Kendall \tau_b)');
title('Ranking stability vs p (IQR band)', 'FontSize', cfg.plot.title_size);
grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);
xlim([min(x)-0.5, max(x)+0.5]);

exportgraphics(fig, [outbase '.png'], 'Resolution', cfg.viz.export_dpi);
if cfg.viz.save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType','vector');
end
close(fig);

end
