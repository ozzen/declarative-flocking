clc
clear;
close all;
addpath('Plotting/');
addpath('m-scripts/');
addpath('m-functions/');
addpath('performance_measures/')

%%

traj_name = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/traj_sim_20.mat';
load(traj_name);

params = traj(1).params;
step_number = 1;
% figure;
% set(gcf, 'Position', [600,100,1200,1200])

dimensions = [4,5];
% trajectory_nums = randi(20, dimensions(1),dimensions(2));
for simNum = 1:prod(dimensions)
    
    subplot(dimensions(1),dimensions(2),simNum);
    
    %% The repeated display code goes here
%     x = traj(simNum).x(step_number,:);
%     y = traj(simNum).y(step_number,:);
%     vx = traj(simNum).vx(step_number,:);
%     vy = traj(simNum).vy(step_number,:);
    
    x = traj(simNum).x;
    y = traj(simNum).y;
    vx = traj(simNum).vx;
    vy = traj(simNum).vy;
%     displayInitState( x, y, vx, vy, 1, [1,0,0]);
    displayTraj( x, y ,vx, vy);

    
    for i = 1:params.n
        for j = i+1:params.n
            edges(i,j) = plot( [x(i), x(j)] , [y(i), y(j)], 'Color',[0 0 0], 'LineWidth', 1.0, 'LineStyle', '-');
            %                 set(edges(i,j), 'Visible', 'off');
            if norm([x(i); y(i)] - [x(j); y(j)]) > params.rs
                set(edges(i,j), 'Visible', 'off');
            else
                set(edges(i,j), 'Visible', 'on');
            end
        end
    end
    
    %%
    axis equal
%     pause(0.5);
    % clf('reset')
end

set(gcf, 'Position', [0 300 1440 1440]);
