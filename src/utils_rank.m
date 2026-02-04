function r = utils_rank(scores)
%UTILS_RANK Rank in descending order with tie handling (average ranks).
%
% Input: scores 1xk
% Output: r 1xk, rank 1 is best (largest score)

scores = scores(:);
k = numel(scores);

[sorted, idx] = sort(scores, 'descend');

r_temp = zeros(k,1);
pos = 1;
while pos <= k
    % find tie block
    val = sorted(pos);
    j = pos;
    while j <= k && sorted(j) == val
        j = j + 1;
    end
    % positions pos..j-1 are tied
    avg_rank = (pos + (j-1))/2;
    r_temp(idx(pos:j-1)) = avg_rank;
    pos = j;
end

r = r_temp(:)';

end
