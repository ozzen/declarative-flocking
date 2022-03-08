function [res] = A_ij(delta_p_ij, params, i, j)
res = zeros(1, params.n * 2);
res(2 * i - 1: 2 * i) = -delta_p_ij;
res(2 * j - 1: 2 * j) =  delta_p_ij;
end

