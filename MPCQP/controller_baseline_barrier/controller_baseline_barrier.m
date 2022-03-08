function [a, fval, exit_flag, history] = controller_baseline_barrier(pos, vel, params, opt)
%% controller_baseline_barrier

%% OutFunction
% Set up shared variables with out_function
history.x = [];
history.fval = [];

opt.OutputFcn = @out_function;

%%
lb = -params.amax * ones(2*params.h,1);
ub = params.amax * ones(2*params.h,1);

u_init = 2 * rand(2 * params.h, 1) - 1;

%% Aeq, beq
 [A, b] = linear_contraints( pos, vel, params);
 
%% fmincon - mpc
% [u, fval, exit_flag] = fmincon(@(u) cost_sum(u, pos, vel, params)...
%     ,u_init,[],[],[],[],lb, ub, @(u) nonlinear_constraints(u, params), opt);

%% fmincon - non-mpc
[u, fval, exit_flag] = fmincon(@(u) cost_lie(u, pos, vel, params)...
    ,u_init,A,b,[],[],lb, ub, @(u) nonlinear_constraints(u, params), opt);

%% fminimax
% [u, fval] = fminimax(@(u) cost_lie(u, pos, vel, params)...
%     ,[rand, rand],[],[],[],[],lb, ub, @(u) nonlinear_constraints(u, params), opt);

%% Solution is processed
a_h = u2acc(u, params);
a = a_h(:,1);
a = trim_vec(a, params.amax);
    
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

