tic
clc
clear;
close all;
source_path = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
addpath([source_path 'Plotting/']);
addpath([source_path 'm-scripts/']);
addpath([source_path 'm-functions/']);
addpath([source_path 'performance_measures/'])
addpath([source_path 'MPCQP/'])

%% Some on/off options
Description = ''; % Added to the output folder name.

selected_agents = [1,8,14];
pair_agents = [2,1]; % For lie-derivative and barrier function plot in simplex. 

skip_1 = 0; %How many steps to skip in video between frames. 

LeaderFollow = 0; 
predator = 0;

PlotEdges = 1; %0:Off, 1:on, 2:selective
PlotEdgeId = selected_agents;

pbf = false;
Obstacles = 0;
pointO = 0;
Target = 0;
target = [60, 50];

Acceleration = 1; %0:Off, 1:on, 2:selective
Acc_scale = 1;
Acc_agent = selected_agents;

Speed = 1; %0:Off, 1:on, 2:selective
speed_agent = selected_agents;

v_scale = 1.0;

simplex_policy = 1;

drawTraj = 0; %1: Incremental, 2: trail
trail_len = 1000;

showVideo = 1;
showText = 1;
FullBox = 0;
BoxL = 2;
showID = 1;
showCircles = 1;

showRegion = 0; % Feasible region of linear barrier based safety constraints. 
agnet_number = 3;

simulation_starting_state = 400; %video starts from

%% Data source, 1: mat file, 2: C code text files
data_source = 1;
if data_source == 1
    traj_num = 1;
    %     experiment = 'bf';
    %     controller = 'cmpc';
    %     traj_name = ['trajs_fossacs2020/' experiment '/' controller '/traj_' controller '_' experiment '_1.mat'];
    %     save_to = ['/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/trajs_fossacs2020/' experiment '/' controller '/'];
    
    traj_path = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/traj/simplex/traj_sim_1 [2020-05-27]_2';
    
    [newStr, matches] = split(traj_path, '/');
    traj_name = [cell2mat(newStr(end)), '.mat']; 
    
    save_to = [traj_path '/Results/'];
    mkdir(save_to);
end

%% Directories Names
directory_setup;
if data_source == 1
    destPath = save_to;
end

%% Load trajectory Data - from matlab code, save files, passes on: x, y, vx, vy, accx, accy, policy
if data_source == 1
    load([traj_path '/' traj_name]);
    x = traj(traj_num).x;
    y = traj(traj_num).y;
    vx = traj(traj_num).vx;
    vy = traj(traj_num).vy;
    accx = traj(traj_num).ax;
    accy = traj(traj_num).ay;
    params = traj(traj_num).params;
%     params.rs = params.rs_bc;
%     mpc_cost = traj.mpc_cost;
%     exit_flag = traj.exit_flag;
%     fitness = traj.fitness; %non-mpc fitness, computed at each time step.
    % params.rad_sensing = params.rs_bc;
    % params.dmin = params.Ds;
    if simplex_policy
        policy = traj(traj_num).policy;
    end
    % bd = traj(traj_num).bd;
    if Obstacles
        [Rect] = readRects(obst_file);
    end
    
    cur_date = datestr(datetime,'(yyyy-mm-dd, hh-MM)');
    dirName = [cur_date ' ' num2str(params.n) '_' num2str(params.steps) '_' num2str(traj_num) '_' Description];
    mkdir([save_to dirName]);
    
%     copyfile(traj_source, [destPath dirName '/'])
elseif data_source == 2
%% Load trajectory Data - from C code
    [params ] = readConf(confFile);

    initFile = sprintf(initPathName, simNum);
    fid = fopen(initFile);
    params.n = str2num(fgets(fid));

    dirName = [num2str(params.n) '_' num2str(params.steps) '_' simNum '_' Description];
    mkdir([destPath dirName]);

    if Obstacles
        [Rect] = readRects(obst_file);
    end

    destination = [destPath dirName '/text_files'];
    saveTextFiles(sourcePath, initPathName, destination, confFile, simNum);
    if Obstacles
        copyfile('rectangles.txt', [destPath dirName '/text_files']);
    end

    [x, y, vx, vy, accx, accy, ~, policy] = read_data(simNum, params.n, sourcePath);
end
%% Save Parameters
% [ret] = saveParams([destPath dirName '/Readme.txt'], params);

%% Set Video options 
videoName = [destPath dirName '/' dirName '-2'];
myVideo = VideoWriter(videoName, 'MPEG-4');
fastForwardRate = 2;
myVideo.FrameRate = min(30,fastForwardRate * (1/params.dt));
myVideo.Quality = 50;
open(myVideo);
skip = skip_1 + 1;

%% Show Graphs
% GraphPlots;

%%
fig_traj = figure;
displayTraj( x, y ,vx, vy);
if Obstacles
    plotRects(obst_file);
end

if pointO
    name = 'pointObstacles25.txt';
    fileID = fopen(name,'r');
    ret = fscanf(fileID,"%f");
    n = ret(1);
    rdata = reshape( ret(2:end) , [2, n] )';
    rdata = rdata(1:n-1,:);
    plot(rdata(:,1),rdata(:,2),'k.', 'MarkerSize', 7);
end

saveas(gcf, [destPath dirName '/' 'trajectory' '.jpg']);
% plot(rdata(:,1),rdata(:,2),'k*', 'MarkerSize', 6);


%%
f1 = figure;
rows = 5;
i = 1;
%% pairwise distance min & max
subplot(rows,2,i);
if predator
    H = PlotDistance(x(:,1:end-1), y(:,1:end-1), params.steps, params );
else
    H = PlotDistance(x, y, params.steps, params );
end
% legend(H, {'Min pairwise distance', 'Max pairwise distance'},'FontSize',17);
axis([0 params.steps 0 5]);
i = i + 1;

%% Acc
subplot(rows,2,i);
plotAcc(accx, accy, params);
i = i + 1;

%% minimum distance - zoomed in veiw
subplot(rows,2,i);
if predator
    min1 = plotMinPWdist(x(:,1:end-1), y(:,1:end-1), params); 
else
    min1 = plotMinPWdist(x, y, params); 
end
set(min1, 'Color', 'r');
plot([1, params.steps], [params.dmin params.dmin], 'g--', 'LineWidth', 1.5);

axis([0 params.steps params.dmin - 0.002 params.dmin + 0.002]);
i = i + 1;

%% Fitness - mpc
% subplot(rows,2,i);
% % plotFileData(sprintf([prefix '/gmpc_%d_fitness.txt'], simNum), 'Fitness Value' ,params.n);
% % i = i + 1;
% if data_source == 1
%     subplot(rows,2,i);
%     plotCurve(mpc_cost, 'Time', 'cost_{mpc}', 1.5, params.n);
%     i = i + 1;
% end

%% Fitness - non mpc
% subplot(rows,2,i);
% % plotFileData(sprintf([prefix '/gmpc_%d_fitness.txt'], simNum), 'Fitness Value' ,params.n);
% % i = i + 1;
% if data_source == 1
%     subplot(rows,2,i);
%     plotCurve(fitness, 'Time', 'Fitness', 1.5, params.n);
%     i = i + 1;
% end

%%  CBF fitness value 
% subplot(rows,2,i)
% prefix = 'experiments/exp1';
% LDfname = sprintf([prefix '/gmpcOBST_%d_LD.txt'], simNum);
% % [res, alpha] = BFbasedFitness( x, y, vx, vy, params, LDfname );
% plotCurve(res, 'Time', 'LD fitness', 1.5, params.n);
% i = i + 1;

%% BF
% subplot(rows,2,i);
% PlotBf(x, y, vx, vy, params);
% % plot(trajMI.h', 'LineWidth', 1.2);
% i = i + 1;

%% Vel
subplot(rows,2,i);
plotVel( vx, vy );
i = i + 1;

%% exit flag for fmincon
% subplot(rows,2,i);
% plotCurve(exit_flag, 'Time', 'exit_flag', 1.5, params.n);
% i = i + 1;

%% Policy
if simplex_policy
    subplot(rows, 2, i);
    plotCurve(policy(:,1)', 'Time', 'Policy', 1.5, params.n);
    ylim([0,3])
    i = i + 1;
end

%% Predator Avoidance
if predator
    subplot(rows,2,i);
    plotPredatorDistance(x, y, params)
    i = i + 1;
end

%% Obstacle Avoidance
if Obstacles
    subplot(rows,2,i);
    plotObstacleDistance(x, y, Rect, params);
    i = i + 1;
end

%% Number of neighbors
% subplot(rows,2,i);
% for rs = [params.rs_bc]
%     num_of_neighbors = zeros(1, params.steps);
%     for index = 1:params.steps
%         q = [x(index,:)', y(index,:)'];
%         pnet = proximityNet(q, rs);
%         num_of_neighbors(index) = numberOfNeighbors(pnet);
%     end
%     plotCurve(num_of_neighbors, 'Time', ['mean num neighbors, rs = ' num2str(rs)], 1.5, params.n);
%     hold on
% end
% i = i + 1;
%% LD
% subplot(rows,2,i);
% plotFileData(LDfname, 'LD' ,params.n);
% hold on 
% plotCurve(-real(alpha), 'Time', 'LD', 1.5, params.n);
% i = i + 1;

%% Learning Rate
% subplot(rows,2,i);
% PlotLearingRate;
% i = i + 1;

%% Learning Rate Log
% subplot(rows,2,i);
% PlotLogLearningRate;
% i = i + 1;


%% Fitness over horizon.
% subplot(rows,2,i);
% figure;
% PlotFitOverH('fitOverHorizon.txt', params)
% i = i + 1;
% % 

%% Barrier funcion between pairs ii and jj
if simplex_policy
    subplot(rows,2,i);
    plotCurve(traj(traj_num).h_ij(pair_agents(1),:,pair_agents(2))...
        , 'Time', ['h_' num2str(pair_agents(1)) '_' num2str(pair_agents(2))], 1.5, params.n);
    i = i + 1;
end

%% Lie-derivative between pairs ii and jj
if simplex_policy
    subplot(rows,2,i);
    plotCurve(traj(traj_num).bd(pair_agents(1),:,pair_agents(2))...
        , 'Time', ['Lie-D_{' num2str(pair_agents(1)) ',' num2str(pair_agents(2)) '}'], 1.5, params.n);
    i = i + 1;
end

%% Acceleration + policy
f2 = figure;
num_plots = numel(selected_agents);
if simplex_policy
    for i = 1:num_plots
        subplot(num_plots, 1, i)
        plotCurve(policy(:,selected_agents(i))', 'Time', '', 1.5, params.n);
        hold on
        plotAcc(accx(2:end,selected_agents(i))', accy(2:end,selected_agents(i))', params)
        title(['Agent ' num2str(selected_agents(i))])
    end
end
%% Save Plots
saveas(f1, [destPath dirName '/' 'combined_plots' '.jpg']);
saveas(f2, [destPath dirName '/' 'accelerations' '.jpg']);

savefig(f1, [destPath dirName '/' 'combined_plots']);
savefig(f2, [destPath dirName '/' 'accelerations']);

%% Show Video
if showVideo
    style = '--';
    traj_color = [0.8 0.8 0.8];
    main2;
end
close(myVideo);
drawnow;



