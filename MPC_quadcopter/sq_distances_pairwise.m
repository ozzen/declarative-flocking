function [res] = sq_distances_pairwise(pos)
%% distances_pairwise - 
% Input:
%   - pos            % 3 x n - state of the flock 
%   - params         % parameters

% Output:
%   - res            % n x n - distances 
% Usama Mehmood - Oct 2019

n = size(pos,2);
res = zeros(n, n);

for i = 1:n
    res(i,:) = sum((pos - repmat(pos(:,i), 1, n)).^2,1);
end

end

