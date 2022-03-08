function [acc] = u2acc(u, params)
%u2acc - convert column vector u to acc matric
% Output:
%   - acc            % n x 2 x h - The sequence of control action over the horizon
%%
acc = reshape(u, [params.n, 2, params.h]);

end

