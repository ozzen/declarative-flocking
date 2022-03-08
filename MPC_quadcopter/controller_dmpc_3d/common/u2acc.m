function [acc] = u2acc(u, params)
%% u2acc - 3D DMPC - convert col vector representation of acc to matrix form. 
% Input:
%    u    - 3.reh x 1 vector
% Output
%    acc    - 3xnxh Matric.
% Usama Mehmood - Oct 2019

acc = reshape(u,[3,params.h]);
end

