function [] = PlotFitOverH(filename, params)

iters = params.maxiters + 1;
range = plotFileData(filename, 'Fitness', params.n);

ys = repmat([range(1) + 0.2* (range(1) - range(2)); range(2) - 0.2* (range(1) - range(2))], 1, params.steps-1);
xs = repmat(iters:iters:iters*(params.steps-1),2,1);
hold on
plot(xs, ys, 'k');

axis([0, (params.steps-1)*(iters), range(2) - 0.2* (range(1) - range(2)), range(1) + 0.2* (range(1) - range(2))])
ax = gca;
ax.XTick = 0:(iters+1)*params.steps/10:(iters+1)*params.steps;
ax.XTickLabel = strsplit(num2str(round(0:params.steps/10:params.steps)));
end

