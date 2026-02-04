function cfg = config()
%CONFIG Configuration for Morris best-p search (deterministic, no noise).
%
% Returns cfg struct with all parameters.

cfg = struct();

% Problem setup
cfg.k = 20;
cfg.P = 4:30;   % 4..30 共 27 个候选



% Step size rule
cfg.delta_rule = "two_step"; % Delta = 2/(p-1)

% Morris budgets
cfg.search.h0 = 120;        % 原来 200，候选多建议降
cfg.search.h_max = 5000;
cfg.search.R0 = 8;          % 原来 10，可略降
cfg.search.eta = 2;      % repeats per p in search stage
cfg.search.min_candidates = 2; % stop if <= this many candidates remain
cfg.search.max_rounds = 10; % safety cap

cfg.final.h = 5000;         % trajectories per run in final evaluation
cfg.final.R = 20;           % repeats per p in final evaluation (finalists only)

% Selection rule thresholds (deterministic => rank should be very stable)
cfg.select.tau_rank = 0.85; % keep p if J_rank >= this
cfg.select.tau_conv = 0.02; % keep p if D_conv <= this (relative change)
cfg.select.eps = 1e-12;

% Separation groups (fixed by design)
cfg.groups.QI = 1:7;
cfg.groups.L  = 8:10;
cfg.groups.W  = 11:20;

% Composite score weights used ONLY for halving (not final decision rule)
cfg.search.score_w = struct('w1',0.6,'w2',0.3,'w3',0.1);

% RNG seeds
cfg.seed.master = 12345;
cfg.seed.beta   = 2026;  % for coefficient generation
cfg.seed.search = 30000; % base seed for search stage runs
cfg.seed.final  = 60000; % base seed for final stage runs

% Output paths
cfg.out.root = pwd;
cfg.out.results_dir = fullfile(cfg.out.root, 'results');
cfg.out.data_dir    = fullfile(cfg.out.results_dir, 'data');



% ---------- Visualization evaluation (uniform across all p) ----------
% 为了画“所有 p 可比”的曲线，建议对所有 p 再做一次统一预算评估（独立于 search/final）
% 默认给一个中等预算：够平滑、又不至于太慢。你要更高质量就调到 5000/20。
cfg.viz = struct();
cfg.viz.enable = true;
cfg.viz.h = 1200;           % 原来 2000，可略降
cfg.viz.R = 8;            % repeats per p for viz dataset
cfg.viz.save_pdf = true; % 同时导出矢量 PDF
cfg.viz.export_dpi = 300;

% ---------- Output dirs ----------
cfg.out.fig_dir = fullfile(cfg.out.results_dir, 'figures');

% ---------- Plot style ----------
cfg.plot = struct();
cfg.plot.font_name = 'Times New Roman';
cfg.plot.font_size = 12;
cfg.plot.title_size = 14;
cfg.plot.line_width = 2.0;
cfg.plot.marker_size = 7;
cfg.plot.grid_alpha = 0.15;
cfg.plot.fig_pos = [100, 100, 1100, 700];  % [x y w h]
cfg.plot.small_pos = [120, 120, 900, 650];

% Group colors (manually chosen, clean + publication-friendly)
cfg.plot.color.QI = [0.20, 0.45, 0.85];  % blue
cfg.plot.color.L  = [0.15, 0.65, 0.35];  % green
cfg.plot.color.W  = [0.55, 0.55, 0.55];  % gray


end
