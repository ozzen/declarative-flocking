function [res] = A_ij(delta_p_ij, params, vi, vj, mode)
% mode 1: normal decentralizded
% mode 2: safety with braking 
if mode == 1
    res = -delta_p_ij;
elseif mode == 2
    if norm(vi) ~= 0
    res = ( norm(vi) * norm(vj) ) / (8 * params.amax^2) * vj +...
        ( dot(vi, vj)*norm(vj) )/ (8 * params.amax^2 * norm(vi)) * vi - ...
        dot(delta_p_ij, vi) / (2 * params.amax * norm(vi)) * vi - ...
        norm(vi) / (2 * params.amax) * delta_p_ij + ...
        params.Ds / params.amax * vi;
    else
    res = [0 0];
    end
end

end

