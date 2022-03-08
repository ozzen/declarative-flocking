function [azi, ele] = vec2azi_ele(n)
%% vec2eze_ele
% This is a test script to compute azimuth and elevation for matlab view(azi, ele)
% values corresponding to a general 3D direction vector n.
% Usama Mehmood June 2020

epsilon = 1e-10;
n = n + epsilon;

azi = 90 + atan360 (n(1), n(2));
ele = atan360(sqrt(n(1)^2 + n(2)^2), n(3));

end

