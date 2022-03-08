function [ res ] = cohesionFitness( x, y )
[Steps, Num] = size(x);
res = zeros(1, Steps);
for s = 1:Num
    for i = 1:Steps
        for j = i+1:Steps
            temp = temp + sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
        end
    end
end
end

