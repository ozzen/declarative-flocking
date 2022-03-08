%% Main script for cmpc experiments
clc
clear
close all

%% dependencies
root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';

addpath([root, 'm-functions']);
addpath([root 'm-functions/init_generators']);

addpath([root, 'MPCQP/common']);

addpath([root, 'MPCQP/controller_dmpc_2d']);
addpath([root, 'MPCQP/controller_dmpc_2d/common']);
addpath([root, 'MPCQP/controller_dmpc_2d/cost_functions']);

%% Local paramters
Runs = 1;
controller = 'dmpc';
experiment = 'bf';

%% Predator params
params.predator = 0;
params.pFactor = 1.40;
params.pred_radius = 6;
params.wp = 380; %350

%% Params:
params.steps = 500;
params.n = 15;
params.h = 3; 
params.dt = 0.1;
params.ct = 0.3;
params.amax = 1.5;
params.vmax = 2;

params.minP = -10;
params.maxP = 10;
params.minV = 0;
params.maxV = 1;
params.Dmax = 100; % neighbor def
params.dmin = 2.0;
params.Ds = params.dmin; % safe dist

params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 0.05; %1e0; %nominal case
% params.rs = sensing_radius(params);
params.rs = 8.4;


params.ws = 30;
params.wc = 1;
params.wo = 125;
params.wt = 1.5;

params.nn = 1;
params.knn = 5;
params.w_m = zeros(1,params.n);


%% Load obstacle file.
rect_file = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/results_neurips/obstacles/1/rectangles.txt';
params.rects = readRects(rect_file);

run controller_dmpc_2d/common/create_mex.m;

%% Optimizer settings
opt = optimoptions('fmincon');
opt.Display = 'off';
% opt.Algorithm = 'active-set';
opt.MaxIterations = 3000;
% opt.StepTolerance = 1e-12;
% Stopping criteria if change in function value is lesser than this:
opt.FunctionTolerance = 1e-3;

%% Load initial states for multiple experiments
% init_file_name = 'istate_15_1000';
% load(['trajs_fossacs2020/src/' init_file_name]);

%%
u = zeros(2*params.n*params.h + 1,1);
indexes = 1:params.n;
acc = zeros(2, params.n);
controller_run = params.ct / params.dt;
zero_vec = zeros(params.n , 1) ;

cpu_time = 0;
function_calls = 0;
tic;

for i = 1:Runs
    %% Load init:
    [posi, veli] = gen_state_random(params);
%         step_number = 1;
%         load('traj_simplex_test.mat');
%         posi = [traj.x(step_number,:) ;traj.y(step_number,:)];
%         veli = [traj.vx(step_number,:) ;traj.vy(step_number,:)];
%         veli = zeros(size(veli));
%     posi = istate(i).pos;
%     veli = istate(i).vel;
    if params.predator 
        posi(:,end) = [30;30];
        veli(:,end) = [0;0];
    end
    %%
    x = zeros([params.steps, params.n]);
    y = zeros([params.steps, params.n]);
    vx = zeros([params.steps, params.n]);
    vy = zeros([params.steps, params.n]);
    ax = zeros([params.steps+1, params.n]);
    ay = zeros([params.steps+1, params.n]);
    cost = zeros(1, params.steps);
    
    x(1,:) = posi(1,:);
    y(1,:) = posi(2,:);
    vx(1,:) = veli(1,:);
    vy(1,:) = veli(2,:);
    
    %% dmpc Controller and Dynamics:
    pos = posi; vel = veli;
    for t = 1:params.steps
        if t == 200
            vel(:,1) = -vel(:,1);
        end
        if mod(t,5) == 0
            e = round(toc, 1);
            disp(['run: ' num2str(i) '/' num2str(Runs) ', step: ' num2str(t) '/' num2str(params.steps) ', Time:' num2str(e) 's, CPU:' num2str(cpu_time)]);
        end
        if mod(t - 1, controller_run) == 0
            if params.predator
                % Controller for predator(nth agent)
                [acc(:,end), dist_to_center] = controller_predator(pos, params);
                if dist_to_center < 10
                    acc(:,end) = [ax(t,end); ay(t,end)];
                end
                for ai = 1:params.n-1 
                    %[posN, velN, neighbors] = get_neighbors_states(ai, pos(:,1:end-1), vel(:,1:end-1), params.rs);
                    [posN, velN, neighbors] = get_k_nearest_neighbors(ai, pos, vel, params.knn);
                    posN = [posN pos(:,end)]; % append predators state
                    velN = [velN vel(:,end)];
                    params.nn = numel(neighbors) + 2;
                    acci = controller_dmpc_2d(posN, velN, params, opt); 
                    acc(:, ai) = acci;             
                end
            else
                for ai = 1:params.n 
%                     [posN, velN, neighbors] = get_neighbors_states(ai, pos, vel, params.rs);
                    [posN, velN, neighbors] = get_k_nearest_neighbors(ai, pos, vel, params.knn);
                    params.nn = numel(neighbors) + 1;
                    
                    if t == 1
                        past_vel = [vx(t, [ai, neighbors]); vy(t, [ai, neighbors])];  
                    else
                        past_vel = [vx(t-1, [ai, neighbors]); vy(t-1, [ai, neighbors])];
                        
                    end
                    
                    delta_vel = sqrt(sum((velN(:,2:end) - past_vel(:,2:end)).^2, 1));
                    
                    if norm(delta_vel) == 0
                        params.w_m(1:numel(neighbors)) = ones(1,numel(neighbors))/numel(neighbors);
                    else
                        params.w_m(1:numel(neighbors)) = delta_vel/sum(delta_vel);                        
                    end
                    
                    start_cpu_time = cputime;
                    acci = controller_dmpc_2d(posN, velN, past_vel, params, opt); 
                    end_cpu_time = cputime;
                    cpu_time = cpu_time + end_cpu_time - start_cpu_time;
                    function_calls = function_calls + 1;
                    
                    acc(:, ai) = acci;             
                end
            end
        else
            acc = [ax(t, :); ay(t, :)]; 
        end
        
        [pos, vel] = dynamics(pos, vel, acc, params);
        x(t+1,:) = pos(1,:);
        y(t+1,:) = pos(2,:);
        vx(t+1,:) = vel(1,:);
        vy(t+1,:) = vel(2,:);
        ax(t+1,:) = acc(1,:);
        ay(t+1,:) = acc(2,:);
    end
    
    %% Store output.
    traj(i).x = x;
    traj(i).y = y;
    traj(i).vx = vx;
    traj(i).vy = vy;
    traj(i).ax = ax(1:params.steps,:);
    traj(i).ay = ay(1:params.steps,:);
    traj(i).params = params;
    
    displayTraj(x,y,vx,vy); title('Minimally Invasive', 'FontSize', 17);
%     set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf, ['Images/trajComparison_' num2str(i) '.jpg']);
%     pause(0.5);
    if mod(i, 10) == 0
    %% Save output
    out_traj_name = ['traj_' controller  '_' experiment '_' num2str(i)];
    out_folder_name = ['trajs_fossacs2020/' experiment '/' controller  '/'];
    save([out_folder_name out_traj_name], 'traj');
    toc;
    end
end

%% Save output
out_traj_name = ['traj_' controller  '_' experiment '_' num2str(Runs)];
out_folder_name = ['trajs_fossacs2020/' experiment '/' controller  '/'];
save([out_folder_name out_traj_name], 'traj');
toc;




