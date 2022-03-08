function [c, ceq] = nonlinear_constraints(u, params) %#codegen
%% contraints - for mpc 
% 1. bound on acceleration magnitude.
% Input:
%    u    - 2 x 1 vector
% Output
%    c    - 1 x 1 vector.
% Usama Mehmood - Oct 2019
%% bound on acceleration magnitude. %TODO for MPC u2acc to be added
    c = sqrt(sum(u.^2,1)) - params.amax;
    ceq = [];
end