function irr = newIrregularity(q, pnet)

irr = 0;
denominator = 0;
[cc,bins] = connectedComponents(pnet);

for comp=1:cc
    filtered_conf = q(bins==comp,:);
    nAgents = size(filtered_conf,1);
    
    if nAgents > 1
        distances = zeros(nAgents,1);
        % for each agent in the CC, find the distance with the closest neighbor
        for i=1:nAgents
            ijs = sqrt(sum((filtered_conf - repmat(filtered_conf(i,:),nAgents, 1)).^2,2));
            ijs(i) = nan;
            distances(i) = min(ijs);
        end
        irr = irr + std(distances,'omitnan');
        denominator = denominator+1;
    end

end

if denominator==0% totally disconnected graph
    irr = 0;
else
    irr = irr/denominator;
end 

end