function [a, fit_val, exit_flag] = controller_dmpc_3d_turn(s, delta_vel, params, opt)
%% controller_dmpc_3d - Run MPC controller and return first action for the entire flock
% Non-linear constraints passed to fmincon as function handle.
% like this:
%     x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
% Input:
%   - s            % 1x12xn - state of the flock
%   - params      d % parameters
%
% Output:
%   - a            % 3 x n - The sequence of control actions
% Usama Mehmood - Oct 2019

a = zeros(3, params.n);
e_delta_vel =  exp(params.eta * delta_vel);

%% find neighbors and arrange inputs to the controller.
for i = 1:params.n
    % Find nearest neighbors
    pos = squeeze(s(1,1:3,:));
    %     vel = squeeze(s(1,4:6,:));
    sq_d = sq_distances_pairwise(pos);
    [~, Index] = sort(sq_d(i,:));
    k_nearest_n = Index(1:params.knn+1);
    
    
    params.w_m(1:numel(k_nearest_n)-1) = e_delta_vel(k_nearest_n(2:end))/sum(e_delta_vel(k_nearest_n(2:end)));
    
    s_local = s(1,:,k_nearest_n);
    %% Optimization
    lb = -params.amax * ones(3*params.h,1);
    ub = params.amax * ones(3*params.h,1);
    
    [u, fit_val, exit_flag] = fmincon(@(u)cost_sum(u, s_local, params),zeros(3*params.h,1),[],[],[],[],lb,ub,@(u)constraints(u, params),opt);
    a_h = u2acc(u, params);
    a(:,i) = a_h(:,1);
    a(:,i) = trim_vec(a(:,i), params.amax);
end

end

