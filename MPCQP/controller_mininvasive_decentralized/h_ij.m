function [res] = h_ij(delta_p_ij, delta_v_ij, vi, vj, params, mode)
% mode 1: normal decentralizded
% mode 2: safety with braking 
if mode == 1
res = sqrt(4 * params.amax * (norm(delta_p_ij) - params.Ds)) + ...
    dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);
elseif mode == 2
    th_i = atan(delta_p_ij / vi);
    th_j = atan(delta_p_ij / vj);
    A1B1_sq = ( (norm(vi)^2 / 4*params.amax) * cos(th_i) - ...
        (norm(vj)^2 / 4*params.amax) * cos(th_j) + ...
        norm(delta_p_ij) )^2 + ...
        ( (norm(vi)^2 / 4*params.amax ) * sin(th_i) - ...
        (norm(vj) / 4*params.amax) * sin(th_j) )^2;
    
    res = A1B1_sq - ...
        (params.Ds + ...
        norm(vi)^2 / 4*params.amax + ...
        norm(vj)^2 / 4*params.amax )^2;  
end

end
