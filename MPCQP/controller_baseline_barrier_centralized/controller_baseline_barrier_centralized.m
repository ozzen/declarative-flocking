function [acc, fval] = controller_baseline_barrier_centralized(pos, vel, params)
%% controller_baseline_barrier_centralized
% Input:
%    pos    - 2 x n vector
%    vel    - 2 x n vector
%    params - struct
% Output
%    acc    - 2 x n vector
%    fval   - 1 x 1 scalar
% Usama Mehmood - Nov 2019

lb = -params.amax * ones(2*params.n, 1);
ub = params.amax * ones(2*params.n, 1);

opt = optimoptions('fmincon');
opt.Display = 'off';

[u, fval] = fmincon(@(u) cost_lie_centralized(u, pos, vel, params)...
,0.3 * ones(2*params.n, 1),[],[],[],[],lb, ub, @(u) nonlinear_constraints_centralized(u, params), opt);

acc = reshape(u, [2, params.n]);
end

