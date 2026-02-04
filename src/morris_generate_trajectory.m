function X = morris_generate_trajectory(k, p, Delta, seed)
%MORRIS_GENERATE_TRAJECTORY Generate a feasible Morris OAT trajectory on p-level grid.
%
% Returns:
%   X: (k+1) x k trajectory points, each in {0,1/(p-1),...,1}
%
% Design:
% - Choose direction s_j in {-1,+1} for each factor
% - Choose starting level so that x_j + s_j*Delta stays within [0,1]
% - Random permutation of factor order for OAT steps

rng(seed);

grid_step = 1/(p-1);

% Sanity: Delta must be multiple of grid_step
ratio = Delta / grid_step;
if abs(ratio - round(ratio)) > 1e-12
    error('Delta=%g is not an integer multiple of grid step=%g for p=%d', Delta, grid_step, p);
end

% Directions
s = randsample([-1, 1], k, true);

% Choose start x0 on grid with feasibility
x0 = zeros(1,k);
max_idx = p-1; % grid indices 0..p-1

% Delta in grid steps (integer)
d_idx = round(ratio);

for j = 1:k
    if s(j) == 1
        % need x_j <= 1-Delta => idx <= max_idx - d_idx
        idx = randi([0, max_idx - d_idx], 1, 1);
    else
        % need x_j >= Delta => idx >= d_idx
        idx = randi([d_idx, max_idx], 1, 1);
    end
    x0(j) = idx * grid_step;
end

% Random order
perm = randperm(k);

X = zeros(k+1, k);
X(1,:) = x0;

x = x0;
for t = 1:k
    j = perm(t);
    x(j) = x(j) + s(j) * Delta;
    % numerical cleanup to grid
    x(j) = round(x(j) / grid_step) * grid_step;
    X(t+1,:) = x;
end

% Final clamp (safety)
X = max(min(X, 1), 0);

end
