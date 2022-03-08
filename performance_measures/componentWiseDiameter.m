function dia = componentWiseDiameter(q,pnet)

dia = 0;

[cc,bins] = connectedComponents(pnet);

for i=1:cc
    filtered_conf = q(bins==i,:);
    distances = pdist(filtered_conf);
    if ~isempty(distances)
        dia = max(dia,max(distances));
    end
end

end