clc
clear 
close all
%%
w = 2;
h = 1;
n = w * (2 * h + 1) - 2 * sum(1:h);
d = 1;
dt = 0.3;
vmax = 1.50;

% pos = zeros(2, n);
vel = zeros(2, n);

jump = sqrt(d^2 - (d/2)^2);
temp = [-0.5*(w-1)*d:d:0.5*(w-1)*d; zeros(1,w)];
pos = temp;

for i = 1:h
    temp = temp(:,1:end - 1);
    temp(1,:) = temp(1,:) + 0.5*d;
    temp(2,:) = temp(2,:) + jump;
    pos = [pos temp];
end

bottom = [pos(1,w+1:end); -pos(2,w+1:end)];
pos = [pos bottom];

plot(pos(1,:), pos(2,:), "r*");
axis equal;
hold on;

for i = 1:n
    point = pos(:,i);
    k = plot(point(1), point(2), 'b*');
    [x, y] = ginput(1);
    delete(k);
    point = [x;y] - point;
    [deg] = atan360(point(1), point(2));
    deg = 60 * floor((deg+30)/60);
    vel(:,i) = [vmax * cosd(deg); vmax * sind(deg)];

end

pos_ = pos + dt*vel;
xD = [pos(1,:) ; pos_(1,:)];
yD = [pos(2,:) ; pos_(2,:)];
l = plot(xD, yD, 'b', 'LineWidth', 1.3); %velocity plot

res = [pos' vel'];
fid = fopen("init_conf_1.txt", 'w');
fprintf(fid, "%d\n", n);
fprintf(fid, "%f %f %f %f\n", res);

save(['state_' num2str(n) '.mat'], 'res');




