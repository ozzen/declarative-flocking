function [ accN, flag, ME, h_allpairs, output, A2, b2] = controller_mininvasive_dec( acc, pos, vel, params)
options = optimoptions('quadprog','Display','off', 'ConstraintTolerance', 1e-1); %ConstraintTolerance: DEFAULT 1E-8
accN = 0; ME = 0;
flag = false;
output = 0;

velSat = vel(:,1) + params.dt * acc;
if norm(velSat) > params.vmax
    velSat = (params.vmax/norm( velSat)) * velSat;
end
acc = (velSat - vel(:,1))/params.dt;

mode = 1;
%Generate Matrices for QP
[H2,f2,A2,b2,lb,ub, h_allpairs] = decentralized_barrier_matrices(pos', vel', acc, params, mode);

%QP
try
    [uNew1, fval, exitflag, output] = quadprog(H2,f2,A2,b2,[],[],lb,ub,[],options);
catch ME
    disp('2d QuadProg has raised exception (Minimally Invasive.)')
    flag = true;
    rethrow(ME);
    return
end

if exitflag == -2
    disp('may day may day');
    flag = true;
    return
end
accN = uNew1; %initialize it.

end

