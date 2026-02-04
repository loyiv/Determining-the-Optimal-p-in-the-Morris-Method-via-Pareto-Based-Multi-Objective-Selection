function plot_mustar_sigma_multiples(cfg, results, outbase)
%PLOT_MUSTAR_SIGMA_MULTIPLES Scatter (mu_star, sigma) per p, with group colors and labels.

P = [results.p];
[~, ord] = sort(P, 'ascend');
results = results(ord);

QI = cfg.groups.QI;
L  = cfg.groups.L;
W  = cfg.groups.W;

for i = 1:numel(results)
    p = results(i).p;
    med_mu_star = results(i).per_run.med_mu_star(:);
    med_sigma   = results(i).per_run.med_sigma(:);

    fig = figure('Visible','off');
    plot_theme(cfg, fig, true);

    hold on;
    scatter(med_mu_star(QI), med_sigma(QI), 55, cfg.plot.color.QI, 'filled');
    scatter(med_mu_star(L),  med_sigma(L),  55, cfg.plot.color.L,  'filled');
    scatter(med_mu_star(W),  med_sigma(W),  55, cfg.plot.color.W,  'filled');

    % factor index labels
    for j = 1:cfg.k
        text(med_mu_star(j), med_sigma(j), sprintf(' %d', j), ...
            'FontSize', 10, 'Color', [0 0 0], 'Clipping', 'on');
    end

    title(sprintf('\\mu^* - \\sigma map (p=%d)', p), 'FontSize', cfg.plot.title_size);
    xlabel('\mu^*'); ylabel('\sigma');
    grid on; set(gca,'GridAlpha',cfg.plot.grid_alpha);

    lg = legend({'QI (1-7)','L (8-10)','W (11-20)'}, 'Location','best');
    lg.Box = 'off';

    exportgraphics(fig, sprintf('%s_p_%d.png', outbase, p), 'Resolution', cfg.viz.export_dpi);
    if cfg.viz.save_pdf
        exportgraphics(fig, sprintf('%s_p_%d.pdf', outbase, p), 'ContentType','vector');
    end
    close(fig);
end
end
