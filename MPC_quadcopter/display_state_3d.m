function [] = display_state_3d(pos, vel, c, params)
v_scale = 1;

pos_ = pos + v_scale * vel;
xD = [pos(1,:) ; pos_(1,:)];
yD = [pos(2,:) ; pos_(2,:)];
zD = [pos(3,:) ; pos_(3,:)];

% c = repmat([1 0 0], params.n, 1);
sze = 300*ones(params.n,1); %Marker Area in points-squared

p = scatter3(pos(1,:), pos(2,:), pos(3,:), sze, c ,'.');
hold on
l = plot3(xD, yD, zD, 'b', 'LineWidth', 1.3); %velocity plot
axis equal

end

