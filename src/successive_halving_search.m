function [finalists, log] = successive_halving_search(cfg, beta)
%SUCCESSIVE_HALVING_SEARCH Successive Halving over p candidates.
%
% Returns:
%   finalists: remaining p values after halving
%   log: struct with per-round metrics and decisions

P_all = cfg.P(:)';
active = P_all;

h = cfg.search.h0;
round_id = 1;

log = struct();
log.rounds = [];

while true
    fprintf('\n[SH] Round %d | Active=%s | h=%d | R0=%d\n', ...
        round_id, mat2str(active), h, cfg.search.R0);

    % Evaluate each active p with current budget
    round_metrics = struct([]);
    round_perrun = struct([]);
    for i = 1:numel(active)
        p = active(i);
        run_seed_base = cfg.seed.search + 10000*round_id + 1000*p;
        [m, pr] = compute_metrics_multirun(cfg, beta, p, h, cfg.search.R0, run_seed_base);
        round_metrics(i).p = p; %#ok<AGROW>
        round_metrics(i).metrics = m; %#ok<AGROW>
        round_perrun(i).per_run = pr; %#ok<AGROW>
    end

    % Compute D_conv for this round among active (sorted ascending)
    [active_sorted, ord] = sort(active, 'ascend');
    % map p -> median mu_star
    med_map = containers.Map('KeyType','double','ValueType','any');
    for i = 1:numel(active)
        p = round_metrics(i).p;
        med_map(p) = round_perrun(i).per_run.med_mu_star;
    end

    Dconv = nan(1, numel(active_sorted));
    for i = 1:numel(active_sorted)-1
        p = active_sorted(i);
        pnext = active_sorted(i+1);
        v = med_map(p);
        vnext = med_map(pnext);
        rel = abs(v - vnext) ./ (abs(vnext) + cfg.select.eps);
        Dconv(i) = max(rel);
    end
    Dconv(end) = NaN;

    % attach D_conv back to each metric
    for i = 1:numel(round_metrics)
        p = round_metrics(i).p;
        pos = find(active_sorted == p, 1);
        round_metrics(i).metrics.D_conv = Dconv(pos);
    end

    % Compute composite score for halving
    Jrank = arrayfun(@(s) s.metrics.J_rank, round_metrics);
    Jsep  = arrayfun(@(s) s.metrics.J_sep,  round_metrics);
    Dcv   = arrayfun(@(s) s.metrics.D_conv, round_metrics);

    % normalize Jsep and Dcv (NaN -> worst)
    Jsep_n = minmax_norm(Jsep);
    Dcv_n  = minmax_norm_nan_worst(Dcv);

    w = cfg.search.score_w;
    S = w.w1 * Jrank + w.w2 * Jsep_n - w.w3 * Dcv_n;

    % sort by score descending
    [~, idx_sort] = sort(S, 'descend');
    active_sorted_by_score = active(idx_sort);

    % log this round
    Rlog = struct();
    Rlog.round = round_id;
    Rlog.h = h;
    Rlog.active = active;
    Rlog.metrics = round_metrics;
    Rlog.score = S;
    Rlog.score_order = active_sorted_by_score;
    log.rounds = [log.rounds; Rlog]; %#ok<AGROW>

    % stopping conditions
    if numel(active) <= cfg.search.min_candidates
        fprintf('[SH] Stop: active candidates <= %d.\n', cfg.search.min_candidates);
        break;
    end
    if h >= cfg.search.h_max
        fprintf('[SH] Stop: h reached h_max=%d.\n', cfg.search.h_max);
        break;
    end
    if round_id >= cfg.search.max_rounds
        fprintf('[SH] Stop: reached max_rounds=%d.\n', cfg.search.max_rounds);
        break;
    end

    % keep top ceil(n/eta)
    n_keep = ceil(numel(active) / cfg.search.eta);
    active = active_sorted_by_score(1:n_keep);
    fprintf('[SH] Keep top %d => %s\n', n_keep, mat2str(active));

    % increase budget
    h = min(h * cfg.search.eta, cfg.search.h_max);
    round_id = round_id + 1;
end

finalists = sort(active, 'ascend');

end

% ---- helpers ----
function y = minmax_norm(x)
xmin = min(x); xmax = max(x);
if abs(xmax - xmin) < 1e-12
    y = zeros(size(x));
else
    y = (x - xmin) / (xmax - xmin);
end
end

function y = minmax_norm_nan_worst(x)
% NaN treated as worst -> set to max before normalization
x2 = x;
mx = max(x(~isnan(x)));
if isempty(mx), mx = 1; end
x2(isnan(x2)) = mx;
y = minmax_norm(x2);
end
