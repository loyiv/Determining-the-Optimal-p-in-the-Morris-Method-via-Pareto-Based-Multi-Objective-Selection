function plot_theme(cfg, fig, isSmall)
%PLOT_THEME Apply consistent publication style.
if nargin < 3, isSmall = false; end

set(fig, 'Color', 'w');

if isSmall
    set(fig, 'Position', cfg.plot.small_pos);
else
    set(fig, 'Position', cfg.plot.fig_pos);
end

set(0, 'DefaultAxesFontName', cfg.plot.font_name);
set(0, 'DefaultTextFontName', cfg.plot.font_name);
set(0, 'DefaultAxesFontSize', cfg.plot.font_size);
set(0, 'DefaultTextFontSize', cfg.plot.font_size);
set(0, 'DefaultLineLineWidth', cfg.plot.line_width);
set(0, 'DefaultLineMarkerSize', cfg.plot.marker_size);
end
