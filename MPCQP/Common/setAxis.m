function [res, limits] = setAxis( pos, params, limits )
    if numel(limits) == 0
        centre = sum(pos,2)/params.n;
        BoxLength = 2.0 * MaxPairwiseDistance(pos(1,:), pos(2,:));
        limits = BoxLength/2*[-1 +1 -1 +1];
        res = centre'*[1 0;1 0;0 1;0 1]' + limits;
    else
        centre = sum(pos,2)/params.n;
        res = centre'*[1 0;1 0;0 1;0 1]' + limits;
    end
end

