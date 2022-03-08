clear;
close all;

addpath('utils')
addpath('performance_measures');

addpath('m-functions');
addpath('m-scripts');
addpath('Plotting')

%% Params
[params ] = readConf( "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/C Code/dmpc.conf" );
params.Tfin = (params.steps - 1) * params.dt;
params.N = 30; % number of agents
params.nd = 2; % space dimension, 2 = 2-D, 3 = 3-D

%% Read Trajs
amax = 1.2;
vmax = 0.75;
beta = 2.5;

nReps = 5;
directory = 'experiments_loop/';
experiment_name = ['exp_' num2string(amax) '_' num2string(vmax) '_' num2string(beta)];

dataFolder = [directory experiment_name];
str1 = 'dmpc_%d_%s.txt';
stateFile = [dataFolder '/' str1];
[ trajs] = makeTrajs( params , nReps, stateFile);
%%
fnams = fieldnames(trajs);
trajs_array = squeeze(struct2cell(trajs));

% 1.'q'              
% 2.'p'              
% 3.'pnet'           
% 4.'connectivity'   
% 5.'v_converg'      
% 6.'irreg'          
% 7.'diameter'       
% 8.'min_pw_dist'    
count = 1;
for i = [3]
    subplot(1,1,count);
    H = printMeasure(fnams{i}, trajs_array, i, params, 'r', '-');
    if i == 3
        plot([0, params.Tfin], [2, 2], 'g--');
        axis([0 params.Tfin 0 8])
        title('Minimum Pairwise distance', 'FontSize', 20);
        legend(H, {'Avg of Min pairwise distance', 'Min of Min pairwise distance'},'FontSize',17);
    end
    count = count + 1;
end

set(gcf, 'Position', [480 300 1440-480 600]);
saveas(gcf, '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Results/MinDist/ws = 0.60.jpg');


[numOfViol, durOfViol] = distanceViolations(trajs, params);
