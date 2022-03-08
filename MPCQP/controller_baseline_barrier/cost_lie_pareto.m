function [cost] = cost_lie_pareto(u, pos, vel, params)
    u = u';
    N = size(pos, 2) - 1;
    cost = [];
    for i = 2:N+1
        delta_p_ij = pos(:, 1) - pos(:, i);
        delta_v_ij = vel(:, 1) - vel(:, i);
        delta_v_bar = dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);

        if delta_v_bar < 0 % Only consider pairs that are collidin

            res = dot(delta_p_ij, u) / norm(delta_p_ij) ...
                - (dot(delta_v_ij, delta_p_ij)^2) / norm(delta_p_ij)^3 * dot(delta_p_ij, vel(:, 1))...
                + dot(delta_v_ij, vel(:, 1)) / norm(delta_p_ij) ...
                + (params.amax * dot(delta_v_ij, delta_p_ij)) / ( norm(delta_p_ij) * sqrt(4 * params.amax * (norm(delta_p_ij) - params.Ds)) );
            
            barrier_function = h_ij(delta_p_ij, delta_v_ij, [], [], params, 1);
            
            cost = [cost, -real(res) / real(barrier_function)];
        end
    end
    
%     distances = sqrt(sum((pos(:,2:end) - pos(:,1)).^2,1));
end

