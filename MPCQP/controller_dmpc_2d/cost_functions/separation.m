function [ret] = separation(pos)
% codegen
%% - separation() - dmpc
n = size(pos, 2);
res = sum((repmat(pos(:,1), 1, n-1) - pos(:,2:end)).^2, 1);
ret = 0;
for i = 1:n-1
    if res(i) ~= 0
        ret = ret + 1/res(i);
    else
        ret = ret + 10000000;
    end
end
ret = ret / (n-1);
end

