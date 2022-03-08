% Zhang 2015 J_v
function J_v = velocityConvergenceNew(p,pnet)

[cc,bins] = connectedComponents(pnet);
J_v = 0;

for comp=1:cc
    filtered_p = p(bins==comp,:);
    compSize = size(filtered_p,1);
    for i=1:compSize
        J_v = J_v +sum((filtered_p(i,:)-mean(filtered_p,1)).^2,2)/compSize;
    end

end
J_v = J_v/cc;
end