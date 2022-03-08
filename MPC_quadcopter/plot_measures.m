clc
clear 
close all

source_path = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/experiments_iccps2020/';

%% Dependencies 
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/performance_measures_3D');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Plotting');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/performance_measures');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/utils');

%% plot
traj_names = {'dnc_quad_bf_noscale.mat', 'dnc_point_bf_noscale.mat', 'cmpc_quad_bf.mat', 'cmpc_point_bf.mat'};
cols = {[1,0,0], [1,0,0], [0,0,1], [0,0,1]};
Marker = {'-', '--', '-', '--'};
field_ids = [22,23,25;...
             14,15,17;...
             21,22,24;...
             13,14,16];

for i = 1:numel(traj_names)
    load([source_path, char(traj_names(i))]);

    fnams = fieldnames(traj(1));
    params = traj(1).params;
    trajs_array = squeeze(struct2cell(traj));
    for m = 1:size(field_ids,2)
        subplot(size(field_ids,2), 1, m)
        H = printMeasure(fnams{field_ids(i,m)}, trajs_array, field_ids(i,m), params, cols{i}, Marker{i});
        hold on
    end
end

legend('Quadcopter/DNC', 'Point model/DNC', 'Quadcopter/CMPC', 'Point model/CMPC', 'FontSize', 17)