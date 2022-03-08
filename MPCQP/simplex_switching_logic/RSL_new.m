function [bool] = RSL_new(pos, vel, params)
%Switch if the Lie Der for all neighbors is greater than x.
N = size(pos, 2);
lb = [-params.amax -params.amax];
ub = [params.amax params.amax];
options = optimoptions(@fmincon, 'Display','off');

bool = true;

for i = 2:N
    delta_p_ij = pos(:, 1) - pos(:, i);
    delta_v_ij = vel(:, 1) - vel(:, i);
    
    [u, fval] = fmincon(@(u) cost_lie_single(u, delta_p_ij, delta_v_ij, params)...
        , [0.3; 0.3],[],[],[],[],lb, ub);
    if -fval < 6 
        bool = false;
        return; 
    end
    %find Maximum value of the lie derivative from all actions.
end

end
