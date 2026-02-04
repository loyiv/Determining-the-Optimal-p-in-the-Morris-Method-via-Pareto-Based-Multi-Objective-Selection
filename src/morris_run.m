function [EE, stats, aux] = morris_run(cfg, beta, p, h, run_seed)
%MORRIS_RUN Run Morris sampling for a given p with h trajectories (deterministic).
%
% Output:
%   EE: h x k matrix
%   stats: struct with mu, mu_star, sigma
%   aux: struct with trajectory seeds (optional)

k = cfg.k;

% Delta rule
Delta = 2/(p-1);

EE = zeros(h, k);

% Generate trajectories and compute EE
for r = 1:h
    traj_seed = run_seed + r; % deterministic per trajectory
    X = morris_generate_trajectory(k, p, Delta, traj_seed);
    Y = model_eval(X, beta); % (k+1)x1
    
    % Each step changes exactly one coordinate => infer which one and compute EE
    for t = 1:k
        x_prev = X(t,:);
        x_next = X(t+1,:);
        diff = x_next - x_prev;
        j = find(abs(diff) > 1e-12);
        if numel(j) ~= 1
            error('Trajectory step does not change exactly one factor at t=%d.', t);
        end
        dy = Y(t+1) - Y(t);
        EE(r, j) = dy / Delta;
    end
end

stats = compute_morris_stats(EE);

aux = struct();
aux.p = p;
aux.h = h;
aux.run_seed = run_seed;
aux.Delta = Delta;

end
