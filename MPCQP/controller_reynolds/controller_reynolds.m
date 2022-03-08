function [a] = controller_reynolds(pos, vel, params, weights)
    w_align = weights.align;
    w_cohes = weights.cohes;
    w_sep = weights.sep;
    
    if size(pos,2) <= 1
        a = [0; 0];
        return
    end
   
    aCohes = mean(pos(:,2:end), 2) - pos(:,1);
    d = sqrt( sum( (pos(:,2:end) - pos(:,1) ).^2 , 1));
    aSep = sum(repmat(1./(d.^2),2,1) .* (pos(:,1) - pos(:,2:end) ),2);
    aAlign =  mean(vel(:,2:end),2) - vel(:,1);
    
    a = w_align*aAlign + w_cohes*aCohes + w_sep*aSep;
%     if norm(a) > params.amax
        a = a * (params.amax / norm(a));
%     end
end