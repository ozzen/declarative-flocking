function [CI_l,CI_u] = confidenceInterval(sampleMeans,sampleStds,n,level)
SEM = sampleStds./sqrt(n);    % Standard Error
levels = sort([level 1-level]);

% lower bound
ts_l = tinv(levels(1),n-1);      % T-Score
CI_l = sampleMeans + ts_l*SEM;    % Confidence Intervals

% upper bound
ts_u = tinv(levels(2),n-1);      % T-Score
CI_u = sampleMeans + ts_u*SEM;    % Confidence Intervals