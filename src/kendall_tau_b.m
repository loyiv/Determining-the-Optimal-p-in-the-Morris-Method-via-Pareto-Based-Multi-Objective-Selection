function tau = kendall_tau_b(rankA, rankB)
%KENDALL_TAU_B Compute Kendall's tau-b correlation between two rankings.
% Handles ties using tau-b correction.
%
% Inputs:
%   rankA, rankB: 1xk rank vectors (smaller rank = better)
% Output:
%   tau in [-1,1]

a = rankA(:);
b = rankB(:);
n = numel(a);

C = 0; D = 0; T_a = 0; T_b = 0;

for i = 1:n-1
    for j = i+1:n
        da = a(i) - a(j);
        db = b(i) - b(j);

        if da == 0 && db == 0
            % tied in both -> ignore in C/D but contributes to tie counts
            T_a = T_a + 1;
            T_b = T_b + 1;
        elseif da == 0
            T_a = T_a + 1;
        elseif db == 0
            T_b = T_b + 1;
        else
            s = sign(da) * sign(db);
            if s > 0
                C = C + 1;
            elseif s < 0
                D = D + 1;
            end
        end
    end
end

den = sqrt((C + D + T_a) * (C + D + T_b));
if den == 0
    tau = 0;
else
    tau = (C - D) / den;
end

end
