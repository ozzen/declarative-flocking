clear;
close all;

addpath('utils')
addpath('performance_measures');

addpath('m-functions');
addpath('m-scripts');
addpath('Plotting')
%% Create and save traj struct for nReps
source_folder = 'experiments_ml/exp29/';

[params ] = readConf( [source_folder 'result_files/gmpc.conf'] );
params.Tfin = (params.steps - 1) * params.dt;
params.N = 30; % number of agents
params.nd = 2; % space dimension, 2 = 2-D, 3 = 3-D

nReps = 25;

str1 = 'gmpc_%d_%s.txt';
stateFile = [source_folder, str1];
% stateFile = 'Shouviks/pred_traj.mat';

[Rect] = readRects([source_folder 'result_files/rectangles.txt']);
target = [60, 50];

% [ traj ] = makeTrajs( params , nReps, stateFile, Rect , target);

res_fldr = [source_folder 'result_files'];
mkdir( res_fldr );

% traj_1 = traj(1:500);
% traj_2 = traj(501:end);
% 
% save([res_fldr '/traj_data_1'], 'traj_1');
% save([res_fldr '/traj_data_2'], 'traj_2');

save([res_fldr '/traj_data'], 'traj');
[ret] = saveParams([res_fldr '/params.txt'], params);

%% Load traj struct and plot results
% source_folder = 'experiments_ml/exp29/';
% res_fldr = [source_folder 'result_files'];
% 
% load([res_fldr '/traj_data']);
% 
% [params ] = readConf( [source_folder 'result_files/gmpc.conf'] );
% params.Tfin = (params.steps - 1) * params.dt;
% params.N = 30; % number of agents
% params.nd = 2; % space dimension, 2 = 2-D, 3 = 3-D

fnams = fieldnames(traj);
trajs_array = squeeze(struct2cell(traj));

% 1.'q'              
% 2.'p' 
% 3. 'a'
% 4. 'obstacle'
% 5.'pnet'           
% 6.'connectivity'   
% 7.'v_converg'      
% 8.'irreg'          
% 9.'diameter'       
% 10.'min_pw_dist'   
% 11. 'min_obstacle_distance'
% 12. 'target'
% 13. 'min_predator_distance'
% 14.  'number_of_pairwise_collisions'
% 15. 'number_of_predator_collisions'
% 16. 'number_of_obstacle_collisions'
% 17. 'avg_number_neighbors'
count = 1;
for i = [6, 7,8, 9, 10, 16]
    subplot(3,2,count);
    H = printMeasure(fnams{i}, trajs_array, i, params, 'r', '-');
    if i == 10
        plot([0, params.Tfin/3], [params.dmin, params.dmin], 'g--');
        axis([0 params.Tfin/3 0 8])
        title('Minimum Pairwise distance', 'FontSize', 20);
        legend(H, {'Avg of Min pairwise distance', 'Min of Min pairwise distance'},'FontSize',17);
    end
    count = count + 1;
end

set(gcf, 'Position', [0 300 1440 600]);

saveas(gcf, [res_fldr '/figures.jpg']);
saveas(gcf, [res_fldr '/figures.fig']);

params.steps = 334;
[numOfViol, durOfViol, res] = distanceViolations(traj, params);

disp(['The number of violations: ' num2str(numOfViol) ' ,duration: ' num2str(durOfViol)]);














