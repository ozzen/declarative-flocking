function [ c , ceq] = myConMI( u, params )
%% 
% Non-liear constraints passed to fmincon as function handle. 
% like this: 
%     x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
% Outputt:
%   - c      % c(x) <= 0 Compute nonlinear inequalities at x.
%   - ceq    % ceq(x) = 0 Compute nonlinear equalities at x.

ceq = 0; 

accelerations = zeros(params.n*params.h, 2);
c = zeros(params.n*params.h, 1);
hor = params.h;
N = params.n;


for k = 0:hor-1
    acc = [u(k*N+1:(k+1)*N) u((hor+k)*N + 1: (hor+k+1)*N)];
    accelerations(k*N+1:(k+1)*N, :) = acc;
end

for m = 1:params.n*params.h
    c(m) = norm(acc(m,:)) - params.amax;
end

end

