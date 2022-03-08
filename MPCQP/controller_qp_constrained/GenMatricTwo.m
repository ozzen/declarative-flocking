function [N ] = GenMatricTwo( n, h )
S = -1 * ones(n) + n * eye(n);
N = zeros(2*n*h+1);
for i = 1:2*h
    N((i-1)*n+1:i*n, (i-1)*n+1:i*n) = S;
end
    N(end, end) = 0;
end

