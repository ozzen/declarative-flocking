function [posN, velN, neighbors] = get_k_nearest_neighbors(ai, pos, vel, knn)
% Find nearest neighbors
%%
distance_to_ai = sqrt(sum((pos - pos(:,ai)).^2, 1));
[~, Index] = sort(distance_to_ai);
k_nearest_n = Index(1:knn+1);

posN = pos(:,k_nearest_n);
velN = vel(:,k_nearest_n);
neighbors = k_nearest_n(2:end);
end

