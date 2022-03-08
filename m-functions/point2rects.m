function [min_distance, intersection_point] = point2rects(point, rects)
% shortest distance between a point and a set of rectangles.
%
% Input:
%   points:             1x2 double
%   rects:              mx1 struct, m rectangles
%
% Output:
%   min_distance:       1x1 double, shortest distance
%   intersection_point: 1x2 double, closest point
m = numel(rects);
min_distance = inf;
intersection_point = [];

for i = 1:m
    [d, intx] = point2rect(point, rects(i));
    if d < min_distance
        min_distance = d;
        intersection_point = intx;
    end
end

end

