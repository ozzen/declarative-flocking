function [s_new] = dynamics_point(s, a ,params)
%% dynamics_point 
% Input:
%   - s      % 1x12xn - state of the flock 
%   - a      % 3xn - The sequence of control action over the horizon 
%   - params % parameters
% Output:
%   - s_new  % next state
% Usama Mehmood - Oct 2019
num_agents = size(s,3);
for i = 1:num_agents
    %update vel
    s(1,4:6,i) = s(1,4:6,i) + params.dt * a(:,i)';
    s(1,4:6,i) = trim_vec(s(1,4:6,i), params.vmax);
    
    %update pos
    s(1,1:3,i) = s(1,1:3,i) + params.dt * s(1,4:6,i);
end
s_new = s;
end

