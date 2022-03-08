function [] = PlotRect(Rect)
    w_vec = 0.5 * Rect.w * [cos(Rect.theta), sin(Rect.theta)];
    beta = pi/2 + Rect.theta;
    l_vec = 0.5 * Rect.l *[cos(beta), sin(beta)];
    mask = [1 1; -1 1; -1 -1; 1 -1];
    corners = repmat([Rect.cx, Rect.cy], 4, 1) + mask*[w_vec; l_vec];
    fill([corners(:,1); corners(1,1)], [corners(:,2); corners(1,2)], [0.5 0.5 0.5], 'LineWidth', 2);
end