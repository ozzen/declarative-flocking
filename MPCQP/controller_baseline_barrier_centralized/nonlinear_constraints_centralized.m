function [c, ceq] = nonlinear_constraints_centralized(u, params)
%% contraints - for mpc 
% 1. bound on acceleration magnitude.
% Input:
%    u    - 3.n.h x 1 vector
% Output
%    c    - n x 1 vector.
% Usama Mehmood - Nov 2019
%% bound on acceleration magnitude.
    acc = reshape(u, [2, params.n]);
    c = sum(acc.^2, 1) - params.amax*params.amax;
    ceq = [];
end