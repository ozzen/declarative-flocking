function [rs] = sensing_radius(params)
%Taken from some paper, Do refer here: []
rs = params.Ds + ...
    (1/ (4 * params.amax)) * ...
    ((4 * params.amax/params.gamma)^(1/3)  + 2*params.vmax)^2;
end

