function [metrics, per_run] = compute_metrics_multirun(cfg, beta, p, h, R, run_seed_base)
%COMPUTE_METRICS_MULTIRUN Run R repeats of Morris at (p,h) and compute metrics.
%
% Outputs:
%   metrics: struct with J_rank, D_conv (placeholder NaN here), J_sep
%   per_run: struct containing per-run stats and ranks

k = cfg.k;

mu_mat      = zeros(R, k);
mu_star_mat = zeros(R, k);
sigma_mat   = zeros(R, k);
rank_mat    = zeros(R, k);

for t = 1:R
    run_seed = run_seed_base + 100*t;
    [EE, stats] = morris_run(cfg, beta, p, h, run_seed);
    mu_mat(t,:) = stats.mu;
    mu_star_mat(t,:) = stats.mu_star;
    sigma_mat(t,:) = stats.sigma;
    rank_mat(t,:) = utils_rank(stats.mu_star);
end

% ---- J_rank: median pairwise Kendall tau-b over ranks
taus = [];
idx = 1;
for a = 1:R-1
    for b = a+1:R
        taus(idx,1) = kendall_tau_b(rank_mat(a,:), rank_mat(b,:)); %#ok<AGROW>
        idx = idx + 1;
    end
end
J_rank = median(taus);

% ---- J_sep: separation in (mu_star, sigma) plane using medians across runs
med_mu_star = median(mu_star_mat, 1);
med_sigma   = median(sigma_mat, 1);

pts = [med_mu_star(:), med_sigma(:)]; % k x 2

QI = cfg.groups.QI(:);
L  = cfg.groups.L(:);
W  = cfg.groups.W(:);

cQI = mean(pts(QI,:), 1);
cL  = mean(pts(L,:), 1);
cW  = mean(pts(W,:), 1);

centers = [cQI; cL; cW];

% group scatter: mean distance to center; then take max
sQI = mean(vecnorm(pts(QI,:) - cQI, 2, 2));
sL  = mean(vecnorm(pts(L,:)  - cL, 2, 2));
sW  = mean(vecnorm(pts(W,:)  - cW, 2, 2));
max_scatter = max([sQI, sL, sW]);

% min inter-center distance
d12 = norm(centers(1,:) - centers(2,:));
d13 = norm(centers(1,:) - centers(3,:));
d23 = norm(centers(2,:) - centers(3,:));
d_min = min([d12, d13, d23]);

J_sep = d_min / (max_scatter + cfg.select.eps);

metrics = struct();
metrics.p = p;
metrics.h = h;
metrics.R = R;
metrics.J_rank = J_rank;
metrics.D_conv = NaN; % filled later when comparing adjacent p
metrics.J_sep = J_sep;

per_run = struct();
per_run.mu = mu_mat;
per_run.mu_star = mu_star_mat;
per_run.sigma = sigma_mat;
per_run.rank = rank_mat;
per_run.med_mu_star = med_mu_star;
per_run.med_sigma = med_sigma;
% median across runs (R x k -> 1 x k)
per_run.med_mu = median(mu_mat, 1);
% optional: mean
% per_run.med_mu = mean(mu_mat, 1);



end
