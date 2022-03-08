clc 
clear all
close all
%%
dir = 'figures/';

xticks = (1:10)';
lineWidth = 1.5;
fig_pos = [400 100 300 200];

dnn_color = [1 0 0];
cmpc_color = [0 0 1];
zero_color = [0 0 0];
colors = [zero_color; cmpc_color; dnn_color];

styles = {'-', '-', '-'};

%% 1 - diameter diff
filename = 'diameter_diff.fig';
figName = 'diff_new_diameter';

openfig([dir filename]);
a = get(gca,'Children');
set(gcf, 'Position', fig_pos)


for i = 1:numel(a)
    set(a(i), 'LineWidth', lineWidth);
    set(a(i), 'Color', colors(i,:));
    set(a(i), 'LineStyle', styles{i});
end

ylim([-0.5, 3.5]);
ylabel('\Delta D', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);

legend('hide');
savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');

%% 2 - vc diff
filename = 'vc_diff.fig';
figName = 'diff_new_vc';

openfig([dir filename]);
a = get(gca,'Children');
set(gcf, 'Position', fig_pos)


for i = 1:numel(a)
    set(a(i), 'LineWidth', lineWidth);
    set(a(i), 'Color', colors(i,:));
    set(a(i), 'LineStyle', styles{i});
end

ylim([-0.5, 1.3]);
ylabel('\Delta VC', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);

legend('hide');
savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');


%% Legend
figName = 'diff_legend';
legend('DNC-DNN', 'CMPC')
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend , 'Orientation', 'horizontal');

savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');
