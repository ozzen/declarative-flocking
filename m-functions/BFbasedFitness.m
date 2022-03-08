function [ res, alpha ] = BFbasedFitness( x, y, vx, vy, params, LDfname )
    dmin = params.dmin;
    maxAcc = params.amax;
    alphaVal = params.alpha;
    LD = ReadFile( LDfname, "%f" );
    LD = LD';
    [Steps, ~] = size(x);
    costBF = zeros(1,Steps);
    for i = 1:Steps
        vel = [vx(i,:) ;vy(i,:)];
        pos = [x(i,:) ;y(i,:)];
        costBF(i) = CompositeBF(pos, vel, dmin, maxAcc );
    end
    alpha = alphaFn(costBF , alphaVal);
    res = 10000*max(-LD - real(alpha), zeros(1, Steps));
end

