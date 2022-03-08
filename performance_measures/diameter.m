function [res] = diameter(pos)

sq_d = sq_distances_pairwise(pos);
res = sqrt(max(max(sq_d)));

end

