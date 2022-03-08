function [a, fit_val, exit_flag, history] = controller_cmpc_3d(s, params, a, opt)
%% controller_cmpc_3d - Run MPC controller and return first action
% Non-liear constraints passed to fmincon as function handle. 
% like this: 
%     x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
% Input:
%   - s            % 1x12xn - state of the flock 
%   - params       % parameters
%   - a            % 3 x n - The sequence of control action over the horizon 

% Output:
%   - a            % 3 x n - The sequence of control action over the horizon 
%   - fit_val      % cost for the optimal sequence of accelerations.
% Usama Mehmood - Oct 2019

%% OutFunction

% Set up shared variables with out_function
history.x = [];
history.fval = [];
opt.OutputFcn = @out_function;

%% Optimization  
lb = -params.amax * ones(3*params.n*params.h,1);
ub = params.amax * ones(3*params.n*params.h,1);

[u, fit_val, exit_flag] = fmincon(@(u)cost_sum(u, s, params),zeros(3*params.n*params.h,1),[],[],[],[],lb,ub,@(u)constraints(u, a, params),opt);
a_h = u2acc(u, params);
a = a_h(:,:,1);
for i = 1:params.n
    a(:,i) = trim_vec(a(:,i), params.amax);
end
    
%% Function to collect the optimizers history. Used in plotting the fval.
function stop = out_function(x,optimValues,state)
    stop = false;
    switch state
        case 'iter'
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x];
        otherwise
    end
end
end

