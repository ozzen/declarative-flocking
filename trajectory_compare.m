clear;
close all;

addpath('utils')
addpath('performance_measures');

addpath('m-functions');
addpath('m-scripts');
addpath('Plotting')
%%
folders = [22, 23, 37, 38, 39];
col = ['r'; 'g'; 'b'; 'k'; 'm']

n = numel(folders);
for i = 1:n
    source_folder = ['experiments_ml/exp', num2str(folders(i)), '/result_files/traj_data.mat'];
    data(i) = load(source_folder);
    
end

[params ] = readConf( ['experiments_ml/exp38/result_files/dmpc.conf'] );
params.Tfin = (params.steps - 1) * params.dt;
params.N = 30; % number of agents
params.nd = 2; % space dimension, 2 = 2-D, 3 = 3-D

%%
fnams = fieldnames(data(1).traj);
measure = 9;

hndles = [];
for i = 1:n
    traj = data(i).traj;
    trajs_array = squeeze(struct2cell(traj));  
    H = printMeasure(fnams{measure}, trajs_array, measure, params, col(i), '-');
    hndles = [hndles H];
end
legend(hndles, {'cmpc', 'dmpc', '5-nn', '7-nn', '9-nn'},'FontSize',17);
ax = gca;
set(ax, 'FontSize', 13);
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




