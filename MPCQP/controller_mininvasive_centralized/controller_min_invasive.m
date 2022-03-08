function [ accN, flag, ME, LD, h_allpairs, output] = controller_min_invasive( acc, pos, vel, params )
    options = optimoptions('quadprog','Display','iter', 'ConstraintTolerance', 1e-1); %ConstraintTolerance: DEFAULT 1E-8
    accN = 0; LD = 0; ME = 0;
    flag = false;
    velSat = vel + params.dt * acc;
    output = 0;
    for j = 1:params.n
        if norm(velSat(:,j)) > params.vmax
            velSat(:,j) = (params.vmax/norm( velSat(:,j) )) * velSat(:,j);
        end
    end
    acc = (velSat - vel)/params.dt;

    %Generate Matrices for QP
    [H2,f2,A2,b2,lb,ub, h_allpairs] = centralized_barrier_matrices(pos', vel', acc, params);
    
    %QP
    try
        [uNew1, fval, exitflag, output] = quadprog(H2,f2,A2,b2,[],[],lb,ub,[],options);
%         if numel(uNew1) == 0
%             options = optimoptions('linprog','Algorithm','dual-simplex');
%             f3 = zeros(size(uNew1));
%             xnew = linprog(f3,A2,b2,[],[],lb,ub,options);
%             disp(xnew);
%         end
    catch ME
        disp('2d QuadProg has messed up(Minimally Invasive.)')
        flag = true;
        A2
        b2
        rethrow(ME);
        return 
    end
    
    if exitflag == -2
        flag = true;
        return
    end
    accN = reshape(uNew1,2, params.n); %initialize it.
%     LD = [A2(1:2:end), A2(2:2:end)] *  [acc(1,:), acc(2,:)]' - b2 + alp; 
end

