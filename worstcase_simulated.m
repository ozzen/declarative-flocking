% clc
% clear all
% close all
%%
f = figure;
% axis([0, 5, -2, 2]);
axis equal
grid off
hold on

posi = [0,0];

v = 1.2;
theta1 = 0;
veli = [v*cosd(theta1), v*sind(theta1)];

plot(posi(1), posi(2), 'g.', 'MarkerSize', 25);
plot([posi(1), posi(1) + 0.1 * veli(1)], [posi(2), posi(2) + 0.1 * veli(2)], 'b-', 'LineWidth', 2.0)

dt = 0.05;
amax = 1;
vmax = 1.50;

for steps = 10:10:100
    temp = zeros(steps,2);
    theta = 45;
    ends_pos = zeros(numel(theta), 2);
    ends_vel = zeros(numel(theta), 2);

    iter = 0;
    for th = theta
        iter = iter + 1;
        acc = [ amax * cosd(th) amax * sind(th) ];
        pos = posi;
        vel = veli;
        for i = 1:steps
            %Update state
            vel = vel + dt * acc;
%             if norm(vel) > vmax
%                 vel = vel * (vmax/norm(vel));
%             end
            pos = pos + dt * vel;
            temp(i,:) = pos;
        end
        ends_pos(iter, :) = temp(end, :);
        ends_vel(iter, :) = vel;
        
        plot([posi(1), temp(:,1)'], [posi(2), temp(:,2)'], 'color', [0.5 0.5 0.5]);
    end
    
    plot(ends_pos(:,1)', ends_pos(:,2)', 'r.', 'MarkerSize', 25);
    plot(ends_pos(:,1)', ends_pos(:,2)', 'b-');

    ends_pos2 = ends_pos + 0.1 * ends_vel;
    plot([ends_pos(:,1) ends_pos2(:,1)]', [ends_pos(:,2) ends_pos2(:,2)]', 'b', 'LineWidth', 1.5);
    plot([ends_pos(:,1) ends_pos(:,1) + 0.90 * amax * cosd(theta)']', [ends_pos(:,2) ends_pos(:,2) + 0.90 * amax * sind(theta)']', 'g', 'LineWidth', 1.5)

    pause(1)
    worstcase_analytic;
end

axis equal

worstcase_analytic;