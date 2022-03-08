function [ acc, fval, LD, alp, Flag ] = controller_QP_constrained( params, pos, vel, u )
%Given Current State i.e. pos, vel, gives Acceleration as output.
    
    Flag = false;
    [H, ~, A, b, Aeq, Beq, alp, lb, ub ] = GenMatQP( params, pos, vel, u);
        try
            [u, fval] = quadprog(H,[],A,b,Aeq,Beq,lb,ub,[],[]);
            
%             [u, fval] = quadprog(H,[],[],[],Aeq,Beq,lb,ub,[],[]);
        catch ME
            disp('QuadProg_Constrained has messed up');
            Flag = true;
            acc = [0; 0] ; fval = 0; LD = 0; alp = 0;
            disp(num2str(A));         
            return;
        end
        if (numel(u) == 0), return; end
        LD = A(1, :) * u - b(1) + alp; 
        acc = u2acc( u, params );
        
%% Back Calculate acceleration.
    [ ~, velSat ] = Dynamics( pos, vel, acc, params );
    acc = (velSat - vel)./params.dt;
end

