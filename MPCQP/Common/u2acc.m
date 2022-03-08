function [ acc ] = u2acc( u, params )
    acc = [u(1:params.n)'; u(params.n*params.h+1:params.n*(params.h+1))'];
end

