function p_star = choose_p_star(cfg, final_results)
%CHOOSE_P_STAR Choose single p_star from finalists using constraints + tie-break.
%
% Rule:
% 1) keep p with J_rank >= tau_rank
% 2) among those, keep p with D_conv <= tau_conv (if D_conv is NaN, treat as failing)
% 3) choose max J_sep
% 4) tie -> smaller p
% If empty after constraints: fallback to max J_rank, then max J_sep, then smallest p.

tau_rank = cfg.select.tau_rank;
tau_conv = cfg.select.tau_conv;

ps = arrayfun(@(s) s.p, final_results);
Jrank = arrayfun(@(s) s.metrics.J_rank, final_results);
Dconv = arrayfun(@(s) s.metrics.D_conv, final_results);
Jsep  = arrayfun(@(s) s.metrics.J_sep,  final_results);

% If D_conv not computed (can happen for single finalist), set D_conv=0
if numel(final_results) == 1
    final_results(1).metrics.D_conv = 0;
    Dconv = 0;
end

ok1 = (Jrank >= tau_rank);
ok2 = (~isnan(Dconv)) & (Dconv <= tau_conv);

cand = find(ok1 & ok2);

if ~isempty(cand)
    % choose max J_sep, tie -> smallest p
    [~, best_rel] = max(Jsep(cand));
    best = cand(best_rel);
    % tie break on p if multiple equal Jsep
    bests = cand(Jsep(cand) == Jsep(best));
    if numel(bests) > 1
        [~, t] = min(ps(bests));
        best = bests(t);
    end
    p_star = ps(best);
    return;
end

% fallback 1: max J_rank
[~, best] = max(Jrank);
bests = find(Jrank == Jrank(best));
if numel(bests) > 1
    % fallback 2: max J_sep among them
    [~, t] = max(Jsep(bests));
    best = bests(t);
    bests2 = bests(Jsep(bests) == Jsep(best));
    if numel(bests2) > 1
        [~, t2] = min(ps(bests2));
        best = bests2(t2);
    end
end

p_star = ps(best);

end
