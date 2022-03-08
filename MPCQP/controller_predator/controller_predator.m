function [a, dist] = controller_predator(pos, params)
%% controller_predator - accelerate towards the centre of the flock
% The nth pos is of the predator.
% Input:
%   - pos          % 2 x n - state of the flock 
%   - params       % parameters

% Output:
%   - a            % 2 x 1 - The  control action of the predator 
% Usama Mehmood - Feb 2020

%% Find the centre
OC = mean(pos(:,1:params.n-1),2);

%% Find the acc
OP = pos(:,params.n);

PC = -OP + OC;
dist = norm(PC);
a = PC * (params.pFactor * params.amax / norm(PC));


end

