function [ ret ] = cost_sum( u, pos, vel, params )
%% cost_sum() - dmpc
% Input:
%   - u            % 2.h x 1
%   - pos          % 2 x n - state of the flock 
%   - vel          % 2 x n - state of the flock 
%   - params       % parameters

% Output:
%   - a            % 2 x n - The sequence of control action over the horizon 
%   - fit_val      % cost for the optimal sequence of accelerations.
% Usama Mehmood - Feb 2020
%%
    pos = pos(:,1:params.nn);
    vel = vel(:,1:params.nn);
    
    acc = u2acc(u, params); %TODO
    control_steps = uint8(params.ct / params.dt);
    res = 0;
    for h = 1:params.h
        a = acc(:,h);
        for k = 1:control_steps
            acceleration = zeros(2, size(pos,2));
            acceleration(:,1) = a;
            [pos, vel] = dynamics(pos, vel, acceleration, params);
        end
        res = res + fitness( pos, params );
    end
%     ret = res;
    ret = res/params.h + sum(u.^2)/params.h;
end

