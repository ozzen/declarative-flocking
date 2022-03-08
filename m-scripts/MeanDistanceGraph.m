% Delaunay
if Num == 2
    avgRadius = sqrt( (x(:,1) - x(:,2)).^2 + (y(:,1) - y(:,2)).^2 );
    figure;
    plot(avgRadius);
    title(['Distance of separation vs time Graph for ' num2str(Num) ' birds. (' num2str(avgRadius(Steps)) ')'],'FontSize', 15);
    xlabel('Time Steps','FontSize', 15);
    ylabel('Distance of Separation','FontSize', 15);
%     set(figure(1), 'Position', get(0, 'Screensize'))
%     saveas(figure(1), ['Results/' dir,'/', dir, '_sep.jpg'])
    return;
end

avgRadius = zeros(Steps, 1);
for index = 1:Steps
    xc = x(index,:);
    yc = y(index,:);
    TRI = delaunay(xc,yc);
    %Calculate average separation.
    total = 0;
    for i = 1:numel(TRI)/3
        temp = TRI(i,:);
        temp1 = [temp(end) temp(1:2)];
        total = total + sum(sqrt((xc(temp) - xc(temp1)).^2 + (yc(temp) - yc(temp1)).^2));
    end
    avgRadius(index) = total/numel(TRI);
end

% figure;
plot(avgRadius, 'LineWidth', 1.5)
hold on 
title(['Average pairwise distance vs time Graph for ' num2str(Num) ' birds.'],'FontSize', 19);
xlabel('Time Steps','FontSize', 18);
ylabel('Distance','FontSize', 18);
ax = gca;
ax.FontSize = 15;
% axis([0 Steps 2.5 3.5])

set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, [destPath dirName '/minDist.jpg'])

%% Delaunay Triangulation Plot. 
figure;
index = Steps;
plot(x(index,:),y(index,:),'r*')
hold on
TRI2 = [TRI TRI(:,1)];
for i = 1:numel(TRI2)/4
    plot(xc(TRI2(i,:)), yc(TRI2(i,:)),'b');
    hold on
end
axis equal
title(['Mean Edge Length = ' num2str(avgRadius(index))],'FontSize', 15)
title(['Triangulation of ' num2str(Num) ' birds. (' num2str(avgRadius(Steps)) ')'],'FontSize', 15);

set(figure(1), 'Position', get(0, 'Screensize'))
% saveas(figure(1), ['Results/' dir,'/', dir, '_tri.jpg'])
