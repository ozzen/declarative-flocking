function [s_new, e, u_dash_prev] = dynamics(s, a, params, e, u_dash_prev)
%% cost_sum 
% Non-liear constraints passed to fmincon as function handle. 
% like this: 
%     x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
% Input:
%   - s      % 1x12xn - state of the flock. 
%   - a      % 3xn - The control action. 
%   - params % parameters
% Output:
%   - cost      % next state
% Usama Mehmood - Oct 2019

num_steps = uint32(params.ct / params.dt);
s_new = zeros(num_steps+1, 12, params.n);
s_new(1,:,:) = s(1 + k*control_jump, :, :);

if params.quad % Plant is a quadrotor
    for i = 1:params.n
        for m = 1:num_steps
            %acceleration to orientation and thrust
            u_dash = acc2thrust(a(:,i), params, u_dash_prev(2:4,i) );
            u_dash_prev(:,i) = u_dash;
            %tracking inner loop control
            [u, e(i)] = controller_inner_pid(s_new(m,:,i), u_dash, e(i));

            [s_new(m+1,:,i), rot_speed(k*control_jump+m,:,i)] = dynamics_quadcopter(s_new(m,:,i), u, params);
        end
    end
else % Plant is point-like
    for m2 = 1:num_steps
        [s_new(m2+1,:,:)] = dynamics_point(s_new(m2,:,:), a ,params);
    end
end


end

