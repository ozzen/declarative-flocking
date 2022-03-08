function [res] = cost_barrier_dec(u, pos, vel, params)
    N = size(pos, 2) - 1;
    res = 0;
    if N == 0
        return;
    end
     
    sum_prev = 0;
    sum_next = 0;
    for i = 2:N+1
        delta_p_ij = pos(:, 1) - pos(:, i);
        delta_v_ij = vel(:, 1) - vel(:, i);
        vi = vel(:, 1);
        vj = vel(:, i); %is correct. dont worry.
        delta_v_bar = dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);
        if delta_v_bar <= 0
            sum_prev = sum_prev + h_ij(delta_p_ij, delta_v_ij, vi, vj, params, 1)^2;
        end
    end
    
    vel(:,1) = vel(:,1) + 20 * params.dt * u;
    if norm(vel(:,1)) > params.vmax
        vel(:,1) = vel(:,1) * (params.vmax / norm(vel(:,1)));
    end
    pos(:,1) = pos(:,1) + 20 * params.dt * vel(:,1);
    
    for i = 2:N+1
        delta_p_ij = pos(:, 1) - pos(:, i);
        delta_v_ij = vel(:, 1) - vel(:, i);
        vi = vel(:, 1);
        vj = vel(:, i); %is correct. dont worry.
        delta_v_bar = dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);
        if delta_v_bar <= 0
            sum_next = sum_next + h_ij(delta_p_ij, delta_v_ij, vi, vj, params, 1)^2;
        end
    end
    
    res = -sum_next + sum_prev;
    
end

