function [c, ceq] = constraints(u, params)
% constraints - non-linear constraints to enforce magnitude of acceleration
%%
acc = u2acc(u, params);
res = sqrt(sum(acc.^2,2)) - params.amax;
c = res(:);
ceq = [];
end

