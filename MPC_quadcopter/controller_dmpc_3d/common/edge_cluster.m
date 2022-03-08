function [agent_ids, acc_turn, rot_axis] = edge_cluster(s, params)
%% edge_cluster - find cluster of agents on the edge
vel = squeeze(s(1,4:6,:));
pos = squeeze(s(1,1:3,:));

avg_vel = mean(vel,2);
centroid = mean(pos,2);

%% Translate the flock to origin
pos_t = pos - centroid;

%% Rotation Matric
rand_dir = [0,0,1]';
dz = avg_vel/norm(avg_vel);

dx = rand_dir - dot(dz, rand_dir) * dz;
dx = dx/norm(dx);

dy = cross(dz, dx);
dy = dy/norm(dy);

R = [dx';dy';dz'];

%% Apply Rotation
pos_r = R*pos_t;
vel_r = R*vel;

%% Farthest from the axis of heading direction
d = vecnorm(pos_r(1:2,:));
[~,I] = max(d);

%% find cluster closest to agent I:

K = convhull(pos_r');
points_on_chull = unique(K(:));

[dd, index] = sort(vecnorm(pos(:,I) - pos));

% agent_ids = zeros(1,params.num_leaders);
% j = 1;
% for i = 1:params.num_leaders
%     while true
%         %         disp(['j:' num2str(j) ', ' num2str(index(j))])
%         %         disp(agent_ids);
%         if any(points_on_chull==index(j))
%             agent_ids(i) = index(j);
%             j = j + 1;
%             break;
%         else
%             j = j + 1;
%         end
%     end
% end

agent_ids = index(1:params.num_leaders);

%% plots
% c = repmat([1 0 0], params.n, 1);
% for ii = agent_ids
%     c(ii,:) = [0,0,0];
% end
% 
% subplot(1,3,1)
% display_state_3d(pos, vel, c, params)
% 
% subplot(1,3,2)
% display_state_3d(pos_t, vel, c, params)
% 
% subplot(1,3,3)
% display_state_3d(pos_r, vel_r, c, params)
% 


%% turn_acc

cluster_centroid = mean(pos(:,agent_ids),2);
d1 = cluster_centroid - centroid;
d1 = d1/norm(d1);
d2 = avg_vel/norm(avg_vel);

dz = cross(d1,d2); % rotation axis
dz = dz/norm(dz);

rand_dir = 2*rand(3,1)-1;

dx = rand_dir - dot(dz, rand_dir) * dz;
dx = dx/norm(dx);

dy = cross(dz, dx);
dy = dy/norm(dy);


R = [dx';dy';dz'];

theta = params.turn_angle;
Rz = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]; %counter-clockwise

acc_turn = R'*Rz*R*avg_vel;
acc_turn = acc_turn/norm(acc_turn);

rot_axis = dz;
% figure;
% display_state_3d(pos, vel, c, params)
% plot3([centroid(1) acc_turn(1)+centroid(1)], [centroid(2) acc_turn(2)+centroid(2)], [centroid(3) acc_turn(3)+centroid(3)], 'k', 'LineWidth', 3)

end

