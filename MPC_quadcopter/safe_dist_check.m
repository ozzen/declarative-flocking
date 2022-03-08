function [res] = safe_dist_check(pos, params)
% Usama Mehmood - Oct 2019

res = true;
pos = squeeze(pos);
for i = 1:params.n
    for j = i+1:params.n
        dist = norm(pos(:,i) - pos(:,j));
        if dist < params.dmin
            res = false;
            return;
        end
    end
end
end

