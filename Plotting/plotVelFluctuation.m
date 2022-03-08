function [] = plotVelFluctuation( vx, vy )
% Velocity fluctuation graph
xLabel = 'Time';
yLabel = 'Fluctuation in vel';
[Steps, Num] = size(vx);
meanVelx =  sum(vx,2)/Num;
meanVely =  sum(vy,2)/Num;
meanSpeed = repmat(sqrt(meanVelx.^2 + meanVely.^2),1, Num);
speed = sqrt(vx.^2 + vy.^2);
speedDev = sqrt((vx - repmat(meanVelx,1,Num)).^2 +  (vy - repmat(meanVely,1,Num)).^2);
res = speedDev ./ (speed + meanSpeed);
plotCurve(res, xLabel, yLabel, 1.3, Num);
axis([0 Steps 0 1.2]);
hold off
end

