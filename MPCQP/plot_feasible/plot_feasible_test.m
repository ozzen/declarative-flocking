%%% Some simple examples to show how to work "plot_feasible.m"

clear;

%% EXAMPLE 1: a plotted region

A = [[2 4];   % 2x + 4y <= 12
     [1 1];   %  x +  y <= 4
     [-1 0];  % non-negativity constraint x>=0
     [0 -1];  % non-negativity constraint y>=0
    ];
b = [12;
     4;
     0;
     0;
    ];
lower_b = [0; 0];
upper_b = [5; 5];
c = [2; 3];


figure(1)
[sorted_vertices, ...
 h_fes, h_bnd, h_fill, h_vert, h_int, h_max, g_labels] = ...
    plot_feasible(A, b, c, lower_b, upper_b, ...
		  'linecolor', 'b', ...
		  'linestyle', '-', ...
		  'filllinestyle', '--', ...
		  'backgroundcolor', [0.6 1 1], ...
		  'linesep', 0.5, ...
		  'plot_vertices', 'ko', ...
		  'label_vertices', 3, ...
		  'label_vertices_size', 14, ...
		  'label_vertices_prec', 0, ...
		  'label_vertices_color', 'b', ...
		  'plot_max', 'r*');
set(h_max, 'markersize', 18)
d = 1.6;
text(d, 4.9, 'maximize', 'fontsize', 18);
text(d, 4.5, '  f = 2x + 3y', 'fontsize', 18, 'FontName', 'Courier');
text(d, 4.1, 'subject to ', 'fontsize', 18);
text(d, 3.7, '  2x + 4y \leq 12', 'fontsize', 18, 'FontName', 'Courier');
text(d, 3.3, '   x +  y \leq 4', 'fontsize', 18, 'FontName', 'Courier');
text(d, 2.9, '        x \geq 0', 'fontsize', 18, 'FontName', 'Courier');
axis square
set(gcf, 'PaperPosition', [0 0 4 4]);
print('-dpng', 'plot_feasible.png');


%% EXAMPLE 2: a plotted region with lots of features turned on

A = [[2 4];   % 2x + 4y <= 12
     [1 1];   %  x +  y <= 4
     [-1 0];  % non-negativity constraint x>=0
     [0 -1];  % non-negativity constraint y>=0
    ];
b = [12;
     4;
     0;
     0;
    ];
lower_b = [0; 0];
upper_b = [7; 7];
c = [2; 3];


figure(2)
[sorted_vertices, ...
 h_fes, h_bnd, h_fill, h_vert, h_int, h_max, g_labels] = ...
    plot_feasible(A, b, c, lower_b, upper_b, ...
		  'linecolor', 'g', ...
		  'linestyle', '-', ...
		  'filllinestyle', '--', ...
		  'backgroundcolor', [1 0.9 0.6], ...
		  'linesep', 0.5, ...
		  'extend_boundaries', ':', ...
		  'plot_intersections', 'rx', ...
		  'plot_vertices', 'ko', ...
		  'label_vertices', 3, ...
		  'label_vertices_size', 10, ...
		  'label_vertices_prec', 1, ...
		  'label_vertices_color', 'b', ...
		  'plot_max', '+');
set(h_max, 'markersize', 15)
legend(h_fill(2), 'iso-objective lines');
axis square



%% EXAMPLE 2: cross-hatching over a region
figure(3)
plot_feasible(A, b, [0;0], lower_b, upper_b, ...
	      'linecolor', 'b', ...
	      'linestyle', '-', ...
	      'linesep', 0.4, ...
	      'lineangle', 80, ...
	      'plot_vertices', 'ko');
% plot_feasible(A, b, c, lower_b, upper_b, ...
% 	      'linecolor', 'b', ...
% 	      'linestyle', '-', ...
% 	      'linesep', 0.4, ...
% 	      'lineangle', -10, ...
% 	      'hold', 1);

axis square


%% EXAMPLE 3: two feasible regions overlapping
figure(4)

A1 = [[2 4];
      [-1 0];  % non-negativity
      [0 -1];  % non-negativity
     ];
b1 = [12;
      0;
      0;
     ];
A2 = [[1 1];
     [-1 0];  % non-negativity
     [0 -1];  % non-negativity
    ];
b2 = [4;
      0;
      0;
     ];
lower_b = [0; 0];
upper_b = [7; 7];
c = [2; 3];

result = plot_feasible(A1, b1, c, lower_b, upper_b, ...
		       'linecolor', 'r', ...
		       'lineangle', 30, ...
		       'linesep', 0.25, ...
		       'hold', 0);

result = plot_feasible(A2, b2, c, lower_b, upper_b, ...
		       'linecolor', 'b', ...
		       'lineangle', -60, ...
		       'linesep', 0.5, ...
		       'hold', 1);

axis square

