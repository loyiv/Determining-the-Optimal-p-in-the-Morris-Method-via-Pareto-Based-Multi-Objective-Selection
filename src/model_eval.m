function y = model_eval(x, beta)
%MODEL_EVAL Deterministic quadratic model evaluation.
%
% Input:
%   x: 1xk or Nxk matrix in [0,1]
% Output:
%   y: scalar or Nx1 vector

if isvector(x)
    x = x(:)'; % 1xk
end

N = size(x,1);
k = size(x,2);

% Linear part
y = beta.beta0 + x * beta.beta1;

% Quadratic pure terms
y = y + sum((x.^2) .* (beta.betaq(:))', 2);

% Interaction terms: 0.5 * x * beta2 * x' with diagonal=0 and symmetric
% Compute efficiently for batch:
% y_int(n) = sum_{i<j} beta2(i,j) x_i x_j
% With symmetric beta2 and diag=0: 0.5 * sum_{i,j} beta2(i,j) x_i x_j
y = y + 0.5 * sum((x * beta.beta2) .* x, 2);

% Ensure column vector
y = reshape(y, [N, 1]);

end
