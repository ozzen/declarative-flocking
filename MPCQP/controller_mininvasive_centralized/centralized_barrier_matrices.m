function [H,f,A,b,lb,ub, h_allpairs] = centralized_barrier_matrices(pos,vel,acc,params)

N = length(pos(:,1));
H = 2*eye(2*N);
f = -2*reshape(acc,2*N,1);

%% Constraints A <= b, Aeq = beq, lb < u < ub
number_of_constraints = nchoosek(N, 2); % number of pairs
A = zeros(number_of_constraints, 2 * N); b = zeros(number_of_constraints, 1); %compute A, b, H, f matries for quadprog
h_allpairs = zeros(number_of_constraints, 1);
%%
count = 1;
for i = 1:N
    for j = i+1:N
        delta_p_ij = pos(i, :) - pos(j, :);
        delta_v_ij = vel(i, :) - vel(j, :);
        
        delta_v_bar = dot(delta_p_ij, delta_v_ij) / norm(delta_p_ij);
        
        if delta_v_bar <= 0
            A(count, :) = A_ij(delta_p_ij, params, i, j);
            b(count) = b_ij(delta_p_ij, delta_v_ij, params);
            h_allpairs(count) = h_ij(delta_p_ij, delta_v_ij, params);
            count = count + 1;
        end
    end
end

lb = -params.amax*ones(2*params.n,1);
ub = params.amax*ones(2*params.n,1);

end
