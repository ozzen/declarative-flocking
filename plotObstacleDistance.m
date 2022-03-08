function [] = plotObstacleDistance(x, y, rects, params)
d = zeros(1, params.steps);
for i = 1:params.steps
    q = [x(i, :)', y(i, :)'];
    [d(i),~,~] = minFlockObstacleDistance(q, rects);
end

plot(1:params.steps, d, 'r-', 'LineWidth', 1.5);
hold on;
plot([0, params.steps], [params.dmin params.dmin], 'g-', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Min obst dist');
title('Obstacle distance')
end

