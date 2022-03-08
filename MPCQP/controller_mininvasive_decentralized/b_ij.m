function [res] = b_ij(delta_p_ij, delta_v_ij, vi, vj, params, mode)
% mode 1: normal decentralizded
% mode 2: safety with braking 
h = h_ij(delta_p_ij, delta_v_ij, vi, vj, params, mode);

if mode == 1
    res = alphaFunction(h, params.gamma) * norm(delta_p_ij) - ...
        (dot(delta_v_ij, delta_p_ij)^2) / norm(delta_p_ij)^2 + ...
        (4 * params.amax * dot(delta_v_ij, delta_p_ij)) / sqrt(4 * params.amax * (norm(delta_p_ij) - params.Ds)) + ...
        norm(delta_v_ij)^2;
    res = 0.5 * real(res);
elseif mode == 2
    res = dot(delta_p_ij, delta_v_ij) + ...
        norm(vi) * dot(delta_v_ij, vi) / (2 * params.amax) + ...
        0.5 * params.gamma * h;
    
end

end