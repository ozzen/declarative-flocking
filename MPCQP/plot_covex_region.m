function [ret] = plot_covex_region(A, b, h, pos)
num_c = numel(b);
if num_c < 3
    return
end

pairs = combnk(1:num_c,2);
solutions = zeros(size(pairs));
for i = 1:size(pairs, 1)
    lin_A = A(pairs(i,:), :);
    lin_b = b(pairs(i,:));
    solutions(i,:) = lin_A\lin_b;
end

solutions(:,1) = solutions(:,1) + pos(1);
solutions(:,2) = solutions(:,2) + pos(2);

if ~ishandle(h)
    ret = fill(solutions(:,1), solutions(:,2), 'r', 'FaceAlpha', 0.2);
else
    set(h, 'XData', solutions(:,1), 'YData', solutions(:,2));
end
end

