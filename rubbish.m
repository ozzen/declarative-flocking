
params.Ds = 1.0;
params.amax = 5.0;
params.vmax = 2.5;

delta_p_ij = params.Ds:0.01:8;
% delta_v_ij = 0:-0.1:-2*params.vmax;
delta_v_ij = params.vmax;

% [X, Y] = meshgrid(delta_p_ij, delta_v_ij);

% 
% res = sqrt(4 * params.amax * (X - params.Ds)) + ...
%     Y ./ X;

res = sqrt(4 * params.amax * (delta_p_ij - params.Ds));% + (2 * params.vmax) ./ delta_p_ij;
% 
% plot(delta_p_ij, res);
% axis equal
% surf(X, Y, res);
plot(delta_p_ij, res);

% xlabel('relative distance');
% ylabel('relative speed');

%%
t = 131;
pos = [x(t,:); y(t,:)];
vel = [vx(t,:); vy(t,:)];
res = wrapper_barrier_function(pos, vel, 1, 2, params)
