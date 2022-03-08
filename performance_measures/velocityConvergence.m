% Zhang 2015 J_v
function J_v = velocityConvergence(p)

N = size(p,1);
J_v = 0;
for i=1:N
    J_v = J_v + sqrt(sum((p(i,:)-mean(p,1)).^2,2));
end

J_v = J_v/N;
end