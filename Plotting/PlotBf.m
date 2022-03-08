function [costBF] = PlotBf(x, y, vx, vy, params)
    [Steps, Num] = size(x);
    costBF = zeros(1,Steps);
    for k = 1:Steps
        vel = [vx(k,:) ;vy(k,:)];
        pos = [x(k,:) ;y(k,:)];
        costBF(k) = CompositeBF(pos, vel, params.dmin, params.amax );
    end
    

    xLabel = 'Time';
    yLabel = 'CBF';
    plotCurve(abs(costBF), xLabel, yLabel, 1.5, Num);
end