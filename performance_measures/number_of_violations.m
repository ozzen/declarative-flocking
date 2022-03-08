function [res] = number_of_violations(x, y, params)
% Inputs:
%       x: steps x agents (double)
%       y: steps x agents (double)
steps = size(x, 1);
n = size(x, 2);

res = 0;

for i = 1:steps
    px = x(i, :);
    py = y(i, :);
    for j = 1:n-1
        
    end
    res = res + violation_step(px, py, params);
end

end

