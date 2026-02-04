function plot_convergence(cfg, results, outbase)
%PLOT_CONVERGENCE Convergence indicator vs p (safe config fields).
% results: struct array with fields p, metrics.D_conv (or filled in main).

P = [results.p];
[~, ord] = sort(P, 'ascend');
results = results(ord);
P = P(ord);

D = nan(1, numel(results));
for i = 1:numel(results)
    D(i) = results(i).metrics.D_conv;
end

% ---- safe style params ----
title_size = get_field(cfg, {'plot','title_size'}, 16);
axis_size  = get_field(cfg, {'plot','axis_size'}, 14);
grid_alpha = get_field(cfg, {'plot','grid_alpha'}, 0.10);
dpi        = get_field(cfg, {'viz','export_dpi'}, 300);
save_pdf   = get_field(cfg, {'viz','save_pdf'}, true);

tau = get_field(cfg, {'select','tau_conv'}, 0.02);

fig = figure('Visible','off');
if exist('plot_theme','file')
    plot_theme(cfg, fig, false);
else
    set(fig, 'Color', 'w');
end

ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

% main line thicker
plot(ax, P, D, '-', 'LineWidth', 2.8);

% reference line thicker
yline(ax, tau, '--', 'LineWidth', 2.2);

% labels
title(ax, 'Convergence indicator vs p', 'FontSize', title_size);
xlabel(ax, 'p', 'FontSize', axis_size);
ylabel(ax, 'D_{conv} (max rel change to next p)', 'FontSize', axis_size);

grid(ax, 'on');
set(ax, 'GridAlpha', grid_alpha);
ax.LineWidth = 1.5;
box(ax, 'on');

% export
exportgraphics(fig, [outbase '.png'], 'Resolution', dpi);
if save_pdf
    exportgraphics(fig, [outbase '.pdf'], 'ContentType', 'vector');
end
close(fig);
end

% -------- helper: safe nested field getter --------
function v = get_field(s, path, default)
v = default;
try
    cur = s;
    for i = 1:numel(path)
        key = path{i};
        if ~isstruct(cur) || ~isfield(cur, key)
            return;
        end
        cur = cur.(key);
    end
    if ~isempty(cur)
        v = cur;
    end
catch
    v = default;
end
end
