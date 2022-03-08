function [H, f, A, b, Aeq, Beq , alp, lb, ub] = GenMatQP( params, pos, vel, u )
    M = GenMatricOne(params.h, pos, vel, params.ct);
    N = GenMatricTwo(params.n, params.h );
%      H = M*N*eye(2*params.n*params.h+1)*N'*M';
    H = M*N*M';

    H = (H + H') / 2;
    f = zeros(2*params.n*params.h + 1, 1);

    [ A, b, alp] = GenAandb( pos, vel, u, params);
    
    Aeq = zeros(1, 2*params.n*params.h + 1);
    Aeq(end) = 1;
    Beq = 1;

    lb = -1.2 * params.amax*ones(2*params.n*params.h+1,1);
    ub = 1.2 * params.amax*ones(2*params.n*params.h+1,1);


end

