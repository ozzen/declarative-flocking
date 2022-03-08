function [] = plot_parameters(matrix, rx, ry, rz)
[ny, nx, nz] = size(matrix);
y = repmat(ry, 1, nx * nz);
x_t = repmat(rx, ny, 1);
x = repmat(x_t(:)', 1, nz);
z_t = repmat(rz, nx * ny, 1);
z = z_t(:)';

number_of_colors = 100;
cls = cool(number_of_colors);
maxim = max(max(max(matrix)));
minim = min(min(min(matrix)));

% c = ones(numel(x), 3) .* repmat(log(matrix(:))/log(maxim), 1, 3);
% test = floor((number_of_colors-1) * log(matrix(:)+1)/log(maxim)) + 1;
test = floor((number_of_colors-1) * matrix(:)/maxim) + 1;

test = min( [test number_of_colors * ones(size(test))], [], 2);
c = cls(test, :);
s = 4000 * ones(size(x));
scatter3(x,y,z,s,c ,'.');

xlabel('vmax');
ylabel('amax');
zlabel('beta');
title('Number of violations');

colormap(cls);
colorbar('Ticks',[0, 1],...
         'TickLabels',{num2str(minim),num2str(maxim)} );

end

