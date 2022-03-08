function [ret] = violation_step(q, params)
n = size(q, 1);
ret = 0;
px = q(:,1)';
py = q(:,2)';
    for j = 1:n
        for k = j+1:n
            dist = norm([px(j) - px(k), py(j) - py(k)]);
            if dist < params.dmin
                ret = ret + 1;
            end
        end
    end
end

