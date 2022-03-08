function [ M ] = GenMatricOne(h, pos, vel, ct)
% input:    h is the size of the prediction horizon.
%           pos is 2 x n matric of position 
%           vel is 2 x n matric of velocities 
n = size(pos, 2);
A = zeros(n*h);
M = zeros(2*n*h+1);
for j = 1:h
    for i = 1:h
        A(n*(j-1)+1:n*j,n*(i-1)+1:n*i) = max(0,i-j+1)*eye(n);
    end
end

A = ct^2 * A; 

trace = zeros(1, n*h);
for m = 1:h
    trace((m-1)*n+1:m*n) = m*ones(1,n);
end
M(1:n*h, 1:n*h) = A;
M(n*h+1:end-1, n*h+1:end-1) = A;
M(end,1:end-1) = [repmat(pos(1,:), 1, h) + ct * trace .* repmat(vel(1,:), 1, h) ...
    , repmat(pos(2,:), 1, h) + ct * trace .* repmat(vel(2,:), 1, h)];
M(end, end) = 0;



end

