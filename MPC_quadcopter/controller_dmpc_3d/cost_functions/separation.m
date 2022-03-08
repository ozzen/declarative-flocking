function [res] = separation(s)
% Usama Mehmood - Oct 2019

pos = squeeze(s(1,1:3,:));
n = size(s, 3);

res = 0;
for i = 1:n
    for j = i+1:n
        res = res + 1/sum((pos(:,i) - pos(:,j)).^2);
    end
end

end

