function [H,f,A,b,lb,ub, h_allpairs] = decentralized_barrier_matrices(pos,vel,acc,params, mode)

N = length(pos(:,1)) - 1;
H = -2*eye(2);
% f = -2*acc';
f = [0 0];

%% Constraints A <= b, Aeq = beq, lb < u < ub
number_of_constraints = N; % number of pairs
A = zeros(number_of_constraints, 2); b = zeros(number_of_constraints, 1); %compute A, b, H, f matries for quadprog
h_allpairs = zeros(number_of_constraints, 1);
%%
count = 1;
for i = 2:N+1
        delta_p_ij = pos(1, :) - pos(i, :);
        delta_v_ij = vel(1, :) - vel(i, :);
        vi = vel(1, :);
        vj = vel(i, :); %is correct. dont worry.
        delta_v_bar = dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);
        
        if delta_v_bar <= 0
            A(count, :) = A_ij(delta_p_ij, params, vi, vj, mode);
            b(count) = b_ij(delta_p_ij, delta_v_ij, vi, vj, params, mode);
            h_allpairs(count) = h_ij(delta_p_ij, delta_v_ij, vi, vj, params, mode);
            count = count + 1;
        end
end

A = A(1:count-1,:);
b = b(1:count-1);

lb = -params.amax*ones(2,1);
ub = params.amax*ones(2,1);

end
