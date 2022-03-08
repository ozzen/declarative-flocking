function [c, ceq] = constraints(u, a, params) %#codegen
%% contraints - for mpc 
% 1. bound on acceleration magnitude.
% 2. bound on maximum deviation in direction of acceleration. 
% Input:
%    u    - 3.n.h x 1 vector
%    a    - 3 x n vector
% Output
%    c    - ~x1 vector.
% Usama Mehmood - Oct 2019
%% bound on acceleration magnitude.
    acc = u2acc(u, params);
    res = sqrt(sum(acc.^2,1)) - params.amax;
    c_1 = res(:);
    ceq = [];
    
%% bound on maximum deviation in direction of acceleration - prediction horizon
%     c_2 = zeros(params.h-1, params.n);
%     for i = 1:params.n
%         for h = 1:params.h-1
%             u = acc(:,i,h);
%             v = acc(:,i,h+1);
%             c_2(h, i) = angle_vectors(u,v) - params.delta_angle;
%         end
%     end
% %% bound on maximum deviation in direction of acceleration - first step
%     c_3 = zeros(params.n, 1);
%     for i = 1:params.n
%         u = a(:,i);
%         v = acc(:,i,1);
%         c_3(i) = angle_vectors(u,v) - params.delta_angle;
%     end  
%     c = cat(1, c_1, c_2(:), c_3);
    c = c_1;

end

