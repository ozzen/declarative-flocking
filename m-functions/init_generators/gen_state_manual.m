function [pos, vel] = gen_state_manual(params, fig)
pos = zeros(2, params.n);
vel = zeros(2, params.n);

%% First agent at origin
% postion 
scatter(pos(1,1), pos(2,1), 1000, 'r', '.');

% max_speed circle
cir = rectangle('Position', [pos(1,1)-params.vmax, pos(2,1)-params.vmax, 2*params.vmax, 2*params.vmax],...
    'LineWidth', 1.0, 'EdgeColor', [0 0 1], 'Curvature', 1);

% Safe distance circle
cir_dmin(1) = rectangle('Position', [pos(1,1)-params.dmin, pos(2,1)-params.dmin, 2*params.dmin, 2*params.dmin],...
    'LineWidth', 1.0, 'EdgeColor', [1 0 0], 'Curvature', 1);

% velocity
[vx, vy] = ginput(1);
x = pos(1,1);
y = pos(2,1);
v = [vx - x, vy - y];
if norm(v) > params.vmax
    v = (params.vmax/norm(v)) * v;
end
vel(:, 1) = v';
plot([x, x + v(1)], [y, y + v(2)], 'b-');

delete(cir);

%% rest of the agents manually placed
figure(fig);
for i = 1:params.n - 1
    [x, y] = ginput(1);
    
    % position update
    scatter(x, y, 1000, 'r', '.');
    pos(:, i+1) = [x; y];
    
    % max_speed circle
    cir_v = rectangle('Position', [x-params.vmax, y-params.vmax, 2*params.vmax, 2*params.vmax],...
        'LineWidth', 1.0, 'EdgeColor', [0 0 1], 'Curvature', 1);
    
    % Safe distance circle
    cir_dmin(i+1) = rectangle('Position', [x-params.dmin, y-params.dmin, 2*params.dmin, 2*params.dmin],...
        'LineWidth', 1.0, 'EdgeColor', [1 0 0], 'Curvature', 1);
    
    % velocity update
    [vx, vy] = ginput(1);
    
    v = [vx - x, vy - y];
    if norm(v) > params.vmax
        v = (params.vmax/norm(v)) * v;
    end
    vel(:, i+1) = v';
    plot([x, x + v(1)], [y, y + v(2)], 'b-');
    delete(cir_v);
end

for i = 1:params.n
    delete(cir_dmin(i));
end
end

