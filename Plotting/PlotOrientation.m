function [] = PlotOrientation( vx, vy )
xLabel = 'Time';
yLabel = 'Orientation';
Num = size(vx,2);
for i = 1:Num
    heading = atan(vy(:,i)./vx(:,i))*180/pi;
    heading = heading + 180*((vy(:,i)<0 & vx(:,i)<0) );
    heading = heading + 180*((vy(:,i)>0 & vx(:,i)<0) ); 
    heading = heading + 360*(heading<0);
    plotCurve(heading, xLabel, yLabel, 1.3, Num);
    hold on
end
hold off 
end

