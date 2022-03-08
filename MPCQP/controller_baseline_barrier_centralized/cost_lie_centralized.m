function [cost] = cost_lie_centralized(u, pos, vel, params)
%% cost_lie_centralized
% Input:
%    u      - 2.n x 1 vector [ax_1 ay_1 ...]
%    pos    - 2 x n vector
%    vel    - 2 x n vector
%    params - struct
% Output
%    cost    - 1 x 1 vector
% Usama Mehmood - Nov 2019

    cost = 0;
    acc = reshape(u, [2, params.n]);
    for i = 1:params.n
        for j = 1:params.n
            if i ~= j
                delta_p_ij = pos(:, i) - pos(:, j);
                delta_v_ij = vel(:, i) - vel(:, j);
                delta_u_ij = acc(:, i) - acc(:, j);
                delta_v_bar = dot(delta_p_ij, delta_v_ij);
                
                if delta_v_bar < 0
                    % Lie-derivatvie from eq. (9) of paper, Safety barrier certificates for heterogenous multi-agent systems. 
                    res = dot(delta_p_ij, delta_u_ij) / norm(delta_p_ij) ...
                    - (dot(delta_v_ij, delta_p_ij)^2) / norm(delta_p_ij)^3  ...
                    + (2 * params.amax * dot(delta_v_ij, delta_p_ij)) / ( norm(delta_p_ij) * sqrt(4 * params.amax * (norm(delta_p_ij) - params.Ds)) ) ...
                    + norm(delta_v_ij)^2 / norm(delta_p_ij);
                    
                    barrier_function  =  h_ij(delta_p_ij, delta_v_ij, [], [], params, 1);
                    cost = cost - real(res) / real(barrier_function);
                end
            end
        end
    end
end