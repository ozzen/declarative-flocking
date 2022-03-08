clc
clear all
close all
%%
dir = 'results_neurips/figures/';

xticks = (1:10)';
lineWidth = 1.5;
fig_pos = [400 100 300 200];

cmpc_color = [0 0 0];
lstm_color = [1 0 0];
dmpc_color = [0 0 1];
dnn_color = [0 0.5 0];
colors = [cmpc_color; lstm_color; dmpc_color; dnn_color ];

styles = {'-', '--', '-.', '--'};

%% 1
filename = 'diam_ca.fig';
figName = 'new_diameter_ca';

openfig([dir filename]);
a = get(gca,'Children');
set(gcf, 'Position', fig_pos)


for i = 1:numel(a)
    set(a(i), 'LineWidth', lineWidth);
    set(a(i), 'Color', colors(i,:));
    set(a(i), 'LineStyle', styles{i});
end

ylim([15, 60]);
ylabel('D', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);

legend('hide');
savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');
%% 2
filename = 'diam_flock.fig';
figName = 'new_diameter_bf';

openfig([dir filename]);
a = get(gca,'Children');
set(gcf, 'Position', fig_pos)


for i = 1:numel(a)
    set(a(i), 'LineWidth', lineWidth);
    set(a(i), 'Color', colors(i,:));
    set(a(i), 'LineStyle', styles{i});
end

ylabel('D', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);

legend('hide');
savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');
%% 3
filename = 'vc_ca.fig';
figName = 'new_vc_ca';

openfig([dir filename]);
a = get(gca,'Children');
set(gcf, 'Position', fig_pos)


for i = 1:numel(a)
    set(a(i), 'LineWidth', lineWidth);
    set(a(i), 'Color', colors(i,:));
    set(a(i), 'LineStyle', styles{i});
end

ylim([0, 1]);
ylabel('VC', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);

legend('hide');
savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');
%% 4
filename = 'vc_flock.fig';
figName = 'new_vc_bf';

openfig([dir filename]);
a = get(gca,'Children');
set(gcf, 'Position', fig_pos)


for i = 1:numel(a)
    set(a(i), 'LineWidth', lineWidth);
    set(a(i), 'Color', colors(i,:));
    set(a(i), 'LineStyle', styles{i});
end

ylim([0, 0.8]);
ylabel('VC', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);

legend('hide');
savefig(gcf, [dir figName '.fig']);
saveas(gcf, [dir figName '.png']);
saveas(gcf, [dir figName '.eps'],'epsc');

legend('show')
%%
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend , 'Orientation', 'horizontal');

