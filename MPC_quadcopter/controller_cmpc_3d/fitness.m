function [res] = fitness(s, params)% Usama Mehmood - Oct 2019
% Usama Mehmood - Oct 2019
res = params.wc * sum_sq_distances(s) + params.ws * separation(s);
end

