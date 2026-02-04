function beta = build_beta(k, seed, groups)
%BUILD_BETA Construct deterministic coefficients for quadratic benchmark model.
%
% Model:
%   Y(x) = beta0 + sum_j beta1(j)*x_j + sum_{i<j} beta2(i,j)*x_i*x_j + sum_j betaq(j)*x_j^2
%
% groups: struct with fields QI, L, W (index arrays)

rng(seed);

beta = struct();
beta.beta0 = 0.0;

% Initialize
beta.beta1 = zeros(k,1);
beta.betaq = zeros(k,1);
beta.beta2 = zeros(k,k); % symmetric, diagonal 0

QI = groups.QI(:);
L  = groups.L(:);
W  = groups.W(:);

% Group L: strong linear main effects, weak quadratic and weak interactions
for jj = L'
    beta.beta1(jj) = unifrnd(8, 12);
    beta.betaq(jj) = unifrnd(-0.2, 0.2);
end

% Group QI: strong quadratic + strong within-group interactions
for jj = QI'
    beta.beta1(jj) = unifrnd(-1, 1);
    beta.betaq(jj) = unifrnd(6, 10);
end

% Group W: weak factors
for jj = W'
    beta.beta1(jj) = unifrnd(-0.1, 0.1);
    beta.betaq(jj) = unifrnd(-0.1, 0.1);
end

% Interactions: within QI strong, cross-group weak
for i = 1:k
    for j = i+1:k
        if ismember(i, QI) && ismember(j, QI)
            v = unifrnd(2, 5);
        elseif ismember(i, W) || ismember(j, W)
            v = unifrnd(-0.1, 0.1);
        else
            % cross QI-L or L-L etc: weak
            v = unifrnd(-0.2, 0.2);
        end
        beta.beta2(i,j) = v;
        beta.beta2(j,i) = v;
    end
end

% Force diagonal to 0
beta.beta2(1:k+1:end) = 0;

end
