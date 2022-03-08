function [A, b  ] = GenConstAcc(params, vel)
    n = params.n;
    h = params.h;
    A = zeros(2 * n * h, 2 * n * h + 1);
    T = zeros(n * h);
    
    for j = 1:h
        for i = 1:j
            T(n*(j-1)+1:n*j,n*(i-1)+1:n*i) = eye(n);
        end
    end
    A(1:n*h, 1:n*h) = T;
    A(n*h+1:end, n*h+1:end-1) = T;
    A = params.ct * A;
    
    b = 0.5*(params.vmax) * ones(2 * n * h, 1) - ...
        [repmat(vel(1,:)', h, 1); repmat(vel(2,:)', h, 1)];
    b = b;
end

