function [res] = b_ij(delta_p_ij, delta_v_ij, params)

h = h_ij(delta_p_ij, delta_v_ij, params);

res = alphaFunction(h, params.gamma) * norm(delta_p_ij) - ...
    (dot(delta_v_ij, delta_p_ij)^2) / norm(delta_p_ij)^2 + ...
    (4 * params.amax * dot(delta_v_ij, delta_p_ij)) / sqrt(4 * params.amax * (norm(delta_p_ij) - params.Ds)) + ...
    norm(delta_v_ij)^2;
res = real(res);
end