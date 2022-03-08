t = steps * dt;
theta = 0:360 ;
r = sqrt(v^2 * t^2 + 0.25 * amax^2 * t^4 + v * amax * cosd(theta) * t^3);

sx = v * cosd(theta) * t + 0.5 * amax * t^2;
sy = v * sind(theta) * t;

alpha = zeros(size(theta));
for i = 1:numel(alpha)
    alpha(i) = atan360(sx(i), sy(i));
end
    
plot(r .* cosd(theta - alpha + theta1), r .* sind(theta - alpha + theta1), 'g-')

axis equal
 