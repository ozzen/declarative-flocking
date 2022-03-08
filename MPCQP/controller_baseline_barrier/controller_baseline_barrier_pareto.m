function [u] = controller_baseline_barrier_pareto(pos, vel, params)
%% controller_baseline_barrier_pareto
    n = size(pos, 2);
    lb = [-params.amax -params.amax];
    ub = [params.amax params.amax];

    opt = optimoptions('fmincon');
    opt.Display = 'off';

    [u1, fval] = fmincon(@(u) cost_lie(u, pos, vel, params)...
    ,[0.3; 0.3],[],[],[],[],lb, ub, @(u) nonlinear_constraints(u, params), opt);
    
    opt = optimoptions('gamultiobj');
    opt.Display = 'off';
    
    [u, fval] = gamultiobj(@(u) cost_lie_pareto(u, pos, vel, params)...
        ,2,[],[],[],[],lb,ub,@(u) nonlinear_constraints(u, params), opt);
    plot(u(:,1), u(:,2), 'r*');
    hold on;
    plot(u1(1), u1(2), 'g*')
    axis equal
    axis([-params.amax +params.amax -params.amax +params.amax]);
    
end

