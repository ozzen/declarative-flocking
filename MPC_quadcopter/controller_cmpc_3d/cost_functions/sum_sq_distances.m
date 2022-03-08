function [res] = sum_sq_distances(s)
% Usama Mehmood - Oct 2019

pos = squeeze(s(1,1:3,:));
n = size(s, 3);

res = 0;
for i = 1:n
    for j = i+1:n
        res = res + sum((pos(:,i) - pos(:,j)).^2);
    end
end


end

