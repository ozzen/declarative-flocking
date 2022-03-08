function [res] = fitness(s, params)% Usama Mehmood - Oct 2019
%% fitness 3D DMPC Usama Mehmood - Oct 2019
if params.turn
    res = params.wc * sum_sq_distances(s) +...
          params.ws * separation(s) +...
          vm(s, params);
else
    res = params.wc * sum_sq_distances(s) + params.ws * separation(s);
end
end

