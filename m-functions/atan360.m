function [res] = atan360(x, y)
%ATAN360 Returns the correct angle from 0 to 360 degrees
% inputs:
%        x: x-coordinate
%        y: y-coordinate
% inputs:
%        ret: angle in degress
% Usama Mehmood - 25th Dec 2018

if x > 0 && y > 0
    deg = atan(y/x);
elseif x <= 0 && y > 0
    deg = pi + atan(y/x);
elseif x <= 0 && y <= 0
    deg = pi + atan(y/x);
elseif x > 0 && y <=0
    deg = 2 * pi + atan(y/x);
end

res = (180/pi) * deg;
end

