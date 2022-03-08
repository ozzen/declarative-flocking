%% main_mpc_quad
% s = zeros(params.steps, 12, params.n);
% accs = zeros(params.steps, 3, params.n);
% s(1,:,:) = s_init;
% 
% control_jump = uint32(params.ct / params.dt);
% control_steps = ceil(params.steps / control_jump);
% 
% lb = -params.amax * ones(3*params.n*params.h,1);
% ub = params.amax * ones(3*params.n*params.h,1);
% a = [zeros(1,params.n); zeros(1,params.n); ones(1,params.n)];
%     
% for i = 1:params.n
%     e(i).prev_e_R = 0;
%     e(i).sum_e_R = 0;
%     e(i).prev_e_w = 0;
%     e(i).sum_e_w = 0;
% end
% u_dash_prev = zeros(4,params.n);
% rot_speed = zeros(params.steps, 4, params.n);
% tic;
% for k = 0:control_steps - 1
%     %% Run MPC controller take first action
%     [u, fval(k+1)] = fmincon(@(u)cost_sum(u, s(1 + k*control_jump,:,:), params),zeros(3*params.n*params.h,1),[],[],[],[],lb,ub,@(u)constraints(u, a, params),opt);
%     a_h = u2acc(u, params);
%     a = a_h(:,:,1);
%     for i = 1:params.n
%         a(:,i) = trim_vec(a(:,i), params.amax);
%     end
%     
%     for c = 0:control_jump-1
%         accs(1+k*control_jump+c, :,:) = reshape(a, [1,3,params.n]);
%     end
%     
%     %% Apply first action to the quadcopters.
% %     [s(1 + k*control_jump:1 + (k+1)*control_jump,:,:)] = dynamics(s(1 + k*control_jump, :, :), a, params);
%     num_steps = uint32(params.ct / params.dt);
%     s_new = zeros(num_steps+1, 12, params.n);
%     s_new(1,:,:) = s(1 + k*control_jump, :, :);
% 
%     if params.quad % Plant is a quadrotor
%         for i = 1:params.n
%             for m = 1:num_steps
%                 
%                 %acceleration to orientation and thrust
%                 u_dash = acc2thrust(a(:,i), params, u_dash_prev(2:4,i) );
%                 u_dash_prev(:,i) = u_dash;
%                 
%                 %tracking inner loop control
%                 [u, e(i)] = controller_inner_pid(s_new(m,:,i), u_dash, e(i));
%                 
%                 [s_new(m+1,:,i), rot_speed(k*control_jump+m,:,i)] = dynamics_quadcopter(s_new(m,:,i), u, params);
%             end
%         end
%     else % Plant is point-like
%         for m2 = 1:num_steps
%             [s_new(m2+1,:,:)] = dynamics_point(s_new(m2,:,:), a ,params);
%         end
%     end
% 
%     s(1 + k*control_jump:1 + (k+1)*control_jump,:,:) = s_new;
%     
%     %% print progress
%     disp(k)
%     
% end
% time_ = toc;
% clc
% clear all
% close all
% %%
% modelfile = 'dnc_model.h5';
% net = importKerasNetwork(modelfile);

input_net = 0.5 * ones(1,36);

res = predict(net,input_net)';
disp(res);







