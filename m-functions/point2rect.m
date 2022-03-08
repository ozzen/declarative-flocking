function [min_distance, intersection_point] = point2rect(point, rect)
% shortest distance between a point and a rectangle.
%
% Input:
%   points:             1x2 double
%   rect:               1x1 struct
%
% Output:
%   min_distance:       1x1 double, shortest distance
%   intersection_point: 1x2 double, closest point

basis_x = [cos(-rect.theta), -sin(-rect.theta)];
basis_y = [sin(-rect.theta), cos(-rect.theta)];

point = point - [rect.cx, rect.cy];
nrm = [dot(point, basis_x), dot(point, basis_y)];

temp = abs(nrm) - [rect.w/2 rect.l/2];

de = max([temp; 0 0], [], 1);

sgn = sign(nrm);

if de(1) ~= 0 && de(2) ~= 0
    intersection_point = [rect.cx, rect.cy] + sgn(1) * rect.w/2 * basis_x +...
        sgn(2) * rect.l/2 * basis_y;
elseif de(1) == 0 && de(2) ~= 0
    intersection_point = [rect.cx, rect.cy] + nrm(1) * basis_x +...
        sgn(2) * rect.l/2 * basis_y;
elseif de(1) ~= 0 && de(2) == 0
    intersection_point = [rect.cx, rect.cy] + sgn(1) * rect.w/2 * basis_x +...
        nrm(2) * basis_y;
else
    intersection_point = [rect.cx, rect.cy] + nrm(1) * basis_x +...
        nrm(2) * basis_y;
end

min_distance = norm(de);
end

