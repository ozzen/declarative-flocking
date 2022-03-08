% note that this version considers std of distances between each agent and
% its neighbors
function irr = paramFreeIrregularity(q,pnet)
N = size(q,1);
% could have as well increased a counter in the 'if ~isempty(q_Ni)'
N_not_isolated = sum(sum(pnet)>0);
irr = 0;
for i=1:N
    q_Ni = q(logical(pnet(i,:)),:);
    if ~isempty(q_Ni)
        distances = sqrt(sum((q_Ni - repmat(q(i,:),size(q_Ni,1), 1)).^2,2));
        irr = irr + std(distances);
    end    
end
irr = irr/N_not_isolated;
end