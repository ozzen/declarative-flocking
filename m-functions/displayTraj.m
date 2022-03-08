function [] = displayTraj( x, y ,vx, vy)
    displayInitState( x, y, vx, vy, 1, [0,1,0])
    plot(x, y, 'LineWidth', 1.0, 'Color', [0.5 0.5 0.5])
    displayInitState( x, y, vx, vy, size(x, 1), [1,0,0])
    axis equal
    set(gcf, 'Position', get(0, 'Screensize'))
end

