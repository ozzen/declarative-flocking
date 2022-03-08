function [] = plotPredatorDistance(x, y, params)
dx = x(:,1:end-1) - x(:,end);
dy = y(:,1:end-1) - y(:,end);
d = sqrt(dx.^2 + dy.^2);

plot(1:params.steps, min(d, [], 2), 'r-', 'LineWidth', 1.5);
hold on
plot(1:params.steps, max(d, [], 2), 'b-', 'LineWidth', 1.5);

plot([0, params.steps],[2, 2], 'g--', 'LineWidth', 1.5);

xlabel('Time/s');
ylabel('Predator distance');
axis([0, params.steps, 0, 5])

end

