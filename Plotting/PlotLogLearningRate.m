range = plotLogFileData('lambda.txt', 'Learning Rate', Num);
ys = repmat([range(1) + 0.2* (range(1) - range(2)); range(2) - 0.2* (range(1) - range(2))], 1, Steps-1);
xs = repmat(iters:iters:iters*(Steps-1),2,1);
hold on
plot(xs, ys, 'k');
axis([0, (Steps-1)*(iters), range(2) - 0.2* (range(1) - range(2)), range(1) + 0.2* (range(1) - range(2))]);