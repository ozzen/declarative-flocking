function [] = plotMaxAcc( vx, vy , Num)
acc = sqrt( ( vx(2:end,:) - vx(1:end-1,:) ).^2 + ( vy(2:end,:) - vy(1:end-1,:) ).^2 ) ;
acc = max(acc, [], 2);
plotCurve(acc, 'Time', 'Max acc', 1.5, Num)
end

