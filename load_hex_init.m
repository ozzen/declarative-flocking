% clc
% clear 
close all
%%
% n = size(res, 1);
n = 4;
load(['state_' num2str(n) '.mat']);
scale = 6.00;
dt = 0.30;
vMax = 1.50;
aMax = 2*1.75;
dmin = 2;

res(:,1:2) = scale * res(:,1:2);
pos = res(:,1:2)';

vel = res(:,3:4)';
% vel = vel * (vMax/2.50);
% res(:,3:4) = vel';

h = plot(pos(1,:), pos(2,:), '.r', 'MarkerSize', 15);
hold on;

pos_ = pos + 3*dt*vel;
xD = [pos(1,:) ; pos_(1,:)];
yD = [pos(2,:) ; pos_(2,:)];
l = plot(xD, yD, 'b', 'LineWidth', 1.3); %velocity plot

for i = 1:n
    posi = pos(:,i)';
    veli = vel(:,i)';
    plotDistanceBrake(posi, veli, vMax, aMax, dt, dmin)
end

axis equal;
%% Plots Lines
% pairs = [1,7];
% for k = 1: size(pairs, 1)
%     xDta = pos(1, pairs(k,:));
%     yDta = pos(2, pairs(k,:));
%     plot(xDta, yDta, 'k');
% end
% %% Move points
% for k = [1,2,4,5,6,7]
%     hold on
%     h2 = plot(pos(1,k), pos(2,k), '.b', 'MarkerSize', 15);
%     [tempx, tempy] = ginput(1);
%     res(k, 1:2) =[tempx, tempy];
%     set(h, 'XData',pos(1, k+1:end), 'YData', pos(2, k+1:end))
%     delete(h2);
%     delete(l(k))
%     plot(tempx, tempy, '.b', 'MarkerSize', 15);
% 
% end

fid = fopen("init_conf_1.txt", 'w');
fprintf(fid, "%d\n", n);

for i = 1:n
    fprintf(fid, "%f %f %f %f\n", res(i,:));
end

%%