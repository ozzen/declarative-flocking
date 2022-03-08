function [res] = average_squared_distance(pos)
% codegen
%% - average_squared_distance() - dmpc
n = size(pos, 2);
res = sum(sum((repmat(pos(:,1), 1, n-1) - pos(:,2:end)).^2, 1), 2);
res = res / (n - 1);
end

