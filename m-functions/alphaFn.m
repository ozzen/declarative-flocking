function [ res ] = alphaFn( costBF , alphaVal)
res = -costBF ./ (costBF - alphaVal);
res = costBF ./ 20;
res = -0.5*ones(size(costBF));
res = 0.25e-05 * costBF.^2;
end

