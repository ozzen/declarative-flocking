function [cost, lie_derivatives, barrier_function] = cost_lie(u, pos, vel, params)
%% cost_lie - The decentralized cost which is the weighted sums of the lie-derivatives of the barrier function in direction of the control input.
% Input:
%   - u            % 2 x 1
%   - pos          % 2 x n - state of the flock 
%   - vel          % 2 x n - state of the flock 
%   - params       % struct of parameters

% Output:
%   - cost                 % 1 x 1 - cost function value 
%   - lie_derivatives      % 1 x n-1 - lie-derivatives 
% Usama Mehmood - 2019

%%
    N = size(pos, 2) - 1;
    ret = zeros(1,N);
    lie_derivatives = zeros(1,N);
    barrier_function = zeros(1,N);
    cost = 0;
    for i = 2:N+1
        delta_p_ij = pos(:, 1) - pos(:, i);
        delta_v_ij = vel(:, 1) - vel(:, i);
        delta_v_bar = dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);

%         if delta_v_bar < 0 % Only consider pairs that are colliding
%             res = dot(delta_p_ij, u) / norm(delta_p_ij) ...
%                 - (dot(delta_v_ij, delta_p_ij)^2) / norm(delta_p_ij)^3 ...
%                 + (2 * params.amax * dot(delta_v_ij, delta_p_ij)) / ( norm(delta_p_ij) * sqrt(4 * params.amax * (norm(delta_p_ij) - params.dmin)) ) ...
%                 + norm(delta_v_ij)^2 / norm(delta_p_ij);

            res = dot(delta_p_ij, u) / norm(delta_p_ij) ...
                - (dot(delta_v_ij, delta_p_ij)) / norm(delta_p_ij)^3 * dot(delta_p_ij, vel(:, 1))...
                + dot(delta_v_ij, vel(:, 1)) / norm(delta_p_ij) ...
                + (params.amax * dot(delta_v_ij, delta_p_ij)) / ( norm(delta_p_ij) * sqrt(4 * params.amax * (norm(delta_p_ij) - params.dmin)) );
            
            lie_derivatives(i-1) = real(res);

            barrier_function(i-1) = h_ij(delta_p_ij, delta_v_ij, [], [], params, 1);
            
%             cost = cost - real(res) / real(barrier_function);
            
            if barrier_function(i-1) ~= 0 %To avoid division by zero.
                ret(i-1) = -lie_derivatives(i-1) / real(barrier_function(i-1));
            else
                ret(i-1) = 100000;
            end
%         end
    end
%     cost = max(ret);
    cost = sum(ret);
end

