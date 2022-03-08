% CreateMap
close all
clear all
%% Read File
res_fldr = 'experiments_ml/exp9/result_files/';
filename = [ 'rectangles.txt'];
%% Plot
[Rect] = readRects(filename);

for i = 1 : size(Rect,2)
    hold on
    PlotRect(Rect(i));
end
axis equal
axis([-20 100 -20 100]);

%%
% point = ginput(1);
% for i = 1:numel(Rect)
%     [d(i), intersection] = point2rect(point, Rect(i));
%     plot(intersection(1), intersection(2), 'r*')
% end
% 
% [~, intersection] = point2rects(point, Rect);
% plot(intersection(1), intersection(2), 'b*')
