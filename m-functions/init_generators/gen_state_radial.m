function [pos, vel] = gen_state_radial(params, radius)
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
% [vx, vy] = ginput(1);
% x = pos(1,1);
% y = pos(2,1);
% v = [vx - x, vy - y];
% if norm(v) > params.vmax
%     v = (params.vmax/norm(v)) * v;
% end
% vel(:, 1) = v';
% plot([x, x + v(1)], [y, y + v(2)], 'b-');
% 
% delete(cir);
vel(:, 1) = [0;0];

%% rest of the agents automatically placed
theta = 0:360/(params.n-1):359;
pos(:,2:end) = [radius * cosd(theta); radius * sind(theta)];
vel(:,2:end) = - pos(:,2:end);

%Display
scatter(pos(1,:), pos(2,:), 1000, 'r', '.');

pos_ = pos + vel;
plot([pos(1,:); pos_(1,:)], [pos(2,:); pos_(2,:)], 'b-');

%% Apply bounds on vel
for i = params.n
    vel(:,2:end) = (params.vmax/norm(vel(:,2:end))) * vel(:,2:end);
end

end

