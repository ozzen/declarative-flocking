function [res] = h_ij(delta_p_ij, delta_v_ij, params)

res = sqrt(4 * params.amax * (norm(delta_p_ij) - params.Ds)) + ...
    dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);

end

