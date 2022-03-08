function [cost] = cost_sum(u, s, params) %#codegen
%% cost_sum - 3D DMPC
% Non-liear constraints passed to fmincon as function handle. 
% like this: 
%     x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
% Input:
%   - s      % 1x12xn - state of the flock 
%   - u      % 3.h x 1 - The sequence of control action over the horizon 
%   - params % parameters
% Output:
%   - cost      % combined cost
% Usama Mehmood - Oct 2019

num_agents = size(s,3);
%% Add zero acc of neighbors
a_i = u2acc(u, params);
a = zeros(3, num_agents, params.h);
a(:,1,:) = a_i;

%%
control_step = params.ct / params.dt;
cost = 0;

e = struct('prev_e_R', repmat({[0;0;0]},1,num_agents),...
            'sum_e_R', repmat({[0;0;0]},1,num_agents),...
            'prev_e_w', repmat({[0;0;0]},1,num_agents),...
            'sum_e_w', repmat({[0;0;0]},1,num_agents));
        
if params.cmpc_prediction_model == 1 %point-model
    for h = 1:params.h
        for k = 1:control_step
            s = dynamics_point(s, a(:,:,h), params);
        end
        cost = cost + fitness(s, params);
    end
elseif params.cmpc_prediction_model == 2 %quadcopter model
    u_dash_prev = zeros(4,num_agents);
    for h = 1:params.h
        for i = 1:num_agents
            for k = 1:control_step
                u_dash = acc2thrust(a(:,i,h), params, u_dash_prev(2:4,i) );
                u_dash_prev(:,i) = u_dash;
                [u_quad, e(i)] = controller_inner_pid(s(1,:,i), u_dash, e(i));
                [s(1,:,i), ~] = dynamics_quadcopter(s(1,:,i), u_quad, params); 
            end
        end
        cost = cost + fitness(s, params);
    end
end

cost = cost/params.h + 0.5*(sum(u.^2)/params.h);

end

