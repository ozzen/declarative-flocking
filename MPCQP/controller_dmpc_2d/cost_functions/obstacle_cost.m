function [res] = obstacle_cost(pos, rects, params)
%% dmpc - obstacle_cost - The cost function for obstacle avoidance.

%%
res = 0;
for j = 1:numel(rects)
    [d, ~] = point2rect(pos(:,1)', rects(j));
    if d ~= 0
        res = res + 1/d^2;
    end
end

count =  numel(rects);
res = res / count;

end

