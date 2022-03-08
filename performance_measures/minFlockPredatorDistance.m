function [minDist, numViolations] = minFlockPredatorDistance(q)

    temp = q(1:end-1, :) - q(end, :);
    dist = sqrt(sum((temp.^2), 2));
    
    numViolations = sum(dist < 1.5);
    minDist = min(dist);
end

