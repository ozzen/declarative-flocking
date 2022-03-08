function [acc] = u2acc(u, params)
%% u2acc - convert col vector representation of acc to matrix form. 
% Input:
%    u    - 3.n.h x 1 vector
% Output
%    acc    - 3xnxh Matric.
% Usama Mehmood - Oct 2019

acc = reshape(u,[3,params.n,params.h]);
end

