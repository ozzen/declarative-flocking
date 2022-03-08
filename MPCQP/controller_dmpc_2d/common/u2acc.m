function [acc] = u2acc(u, params)
%% u2acc - dmpc - convert column vector u to acc matric
% Output:
%   - acc            % 2 x h - The sequence of control action over the horizon
%%
acc = reshape(u, [params.h, 2])';

end

