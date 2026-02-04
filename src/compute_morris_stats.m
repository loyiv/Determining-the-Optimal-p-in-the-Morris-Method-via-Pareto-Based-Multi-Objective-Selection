function stats = compute_morris_stats(EE)
%COMPUTE_MORRIS_STATS Compute Morris indices from EE samples.
%
% Input:
%   EE: h x k matrix of elementary effects
% Output:
%   stats.mu, stats.mu_star, stats.sigma (1 x k)

h = size(EE,1);

stats = struct();
stats.mu = mean(EE, 1);
stats.mu_star = mean(abs(EE), 1);

if h <= 1
    stats.sigma = zeros(1, size(EE,2));
else
    stats.sigma = std(EE, 0, 1); % sample std (normalization by h-1)
end

end
