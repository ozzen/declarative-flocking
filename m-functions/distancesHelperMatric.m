function [ M ] = distancesHelperMatric( n )
M = -ones(n,n) + diag(2*ones(n,1));
end

