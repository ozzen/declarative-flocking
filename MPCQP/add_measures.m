clc
clear all
close

%% Dependencies 
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Plotting');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-functions');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/performance_measures');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/utils');

%% Load traj
% obstacle_field_id = 4;
% controller = 'cmpc';
% experiment = 'bf';
% source_path = ['/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/trajs_fossacs2020/' experiment '/' controller '/'];
% traj_name = ['traj_' controller '_oa_25_' num2str(obstacle_field_id) '.mat'];
% traj_name = ['traj_' controller '_' experiment '_100.mat'];

source_path = '';
traj_name = 'traj_sim_20 [2020-05-16].mat';
load([source_path, traj_name]);
params = traj(1).params;
params.predator = 0;

%% Load rects
% rect_file = ['/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/results_neurips/obstacles/' num2str(obstacle_field_id) '/rectangles.txt'];
% rects = readRects(rect_file);

%% Add measures and update traj
nReps = numel(traj);
% 
for i = 1:nReps
    traj(i).diameter = zeros(params.steps, 1);
    traj(i).connectivity = zeros(params.steps, 1);
    traj(i).velocity_convergence = zeros(params.steps, 1);
    traj(i).min_pw = zeros(params.steps, 1);
%     traj(i).pnet = zeros(params.steps, params.n, params.n);
%       traj(i).obstacle = zeros(params.steps,params.n,2);
%       traj(i).number_of_obstacle_collisions = zeros(params.steps, 1);
%       traj(i).min_obstacle_distance = zeros(params.steps, 1);
%         traj(i).number_of_predator_collisions = zeros(params.steps, 1);
%         traj(i).min_predator_distance = zeros(params.steps, 1);



    for j = 1:traj(1).params.steps
        pos = [traj(i).x(j,:); traj(i).y(j,:)];
        vel = [traj(i).vx(j,:); traj(i).vy(j,:)];
        
        pnet = proximityNet(pos', params.rs);
        traj(i).pnet(j,:,:) = pnet;
        traj(i).connectivity(j) = connectedComponents(pnet);
        if params.predator
            traj(i).min_pw(j) = MinPairwiseDistance(traj(i).x(j,1:end-1), traj(i).y(j,1:end-1));
            traj(i).diameter(j) = diameter(pos(:,1:end-1));
            traj(i).velocity_convergence(j) = velocityConvergence(vel(:,1:end-1));
        else
            traj(i).min_pw(j) = MinPairwiseDistance(traj(i).x(j,:), traj(i).y(j,:));
            traj(i).diameter(j) = diameter(pos);
            traj(i).velocity_convergence(j) = velocityConvergence(vel');
        end
%         [traj(i).min_obstacle_distance(j),...
%          traj(i).obstacle(j, :, :),...
%          traj(i).number_of_obstacle_collisions(j)] = minFlockObstacleDistance(pos', rects);
%         [traj(i).min_predator_distance(j),...
%          traj(i).number_of_predator_collisions(j)] = minFlockPredatorDistance(pos');
%   
    end
    
    if mod(i,5) == 0
        disp(num2str(i));
    end
end
save([source_path, traj_name], 'traj', '-v7.3');

%% Plot
fne = [];
fnams = fieldnames(traj(1));
params = traj(1).params;
trajs_array = squeeze(struct2cell(traj));

fns = {'min_pw', 'diameter', 'velocity_convergence'};


for i = 1:numel(fns)
    fne(i) = fieldName2id(fnams, fns{i});
end

for m = 1:numel(fne)
    subplot(numel(fne), 1, m)
    H = printMeasure(fnams{fne(m)}, trajs_array, fne(m), params, 'r', '-');
    if m == 1
        plot([0, params.steps*params.dt], [params.dmin, params.dmin], 'g--');
    end
%     ylim(gca, [0, 5]);
%     ylabel('distance');
%     legend(H, {'mean of the distance', 'minimum of the distance'})
%     hold on
end

%% Minimum distance constraint violations
[numOfViol, durOfViol, res] = distanceViolations(traj, params);

%% Save fig
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, [source_path  'graphs.jpg']);



% collisions = sum([traj.number_of_predator_collisions],1);

