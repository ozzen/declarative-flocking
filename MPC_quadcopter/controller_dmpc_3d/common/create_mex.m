%% code gen cost_sum(u, s, params)

params_type = coder.typeof(params);
u_type = coder.typeof(zeros(3*params.h,1));
s_type = coder.typeof(zeros(1,12,params.knn+1));

codegen cost_sum.m -args {u_type, s_type, params_type} -o cost_sum

%% code gen angle_vectors(u,v)

% vec_type = coder.typeof([0;0;0]);
% 
% codegen angle_vectors.m -args {vec_type, vec_type} -o angle_vectors

%% code gen constraints(u, a, params)

u_type = coder.typeof(zeros(3*params.h,1));
params_type = coder.typeof(params);

codegen constraints.m -args {u_type, params_type} -o constraints



