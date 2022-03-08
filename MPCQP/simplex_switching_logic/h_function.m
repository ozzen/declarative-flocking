function [ret, delta_v] = h_function(pos_rel, vel_rel, params)

p_rel_n = pos_rel ./ repmat(sqrt(sum(pos_rel.^2, 1)),2 , 1);

delta_v =  - sum(p_rel_n .* vel_rel, 1);

ret = 2 * sqrt(params.amax * (sqrt(sum(pos_rel.^2, 1)) - params.Ds)) -...
    delta_v;
end

