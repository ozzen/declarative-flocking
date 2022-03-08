function C = mouseMove (object, eventdata, h, point, scale)
C = get (gca, 'CurrentPoint');

vec = [C(1,1) - point(1), C(1,2) - point(2)];
magnitude = norm(vec);
if magnitude > scale
    vec = (scale/norm(vec)) * vec;
end
xD = [point(1), point(1) + vec(1)];
yD = [point(2), point(2) + vec(2)];

set(h, 'XData', xD, 'YData', yD);
title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
