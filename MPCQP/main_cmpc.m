%% Main script for cmpc experiments
clc
clear;
close all;

%% dependencies
root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';

addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-functions');
addpath([root 'm-functions/init_generators']);

addpath([root, 'MPCQP/common']);

addpath([root, 'MPCQP/controller_cmpc_2d']);
addpath([root, 'MPCQP/controller_cmpc_2d/common']);
addpath([root, 'MPCQP/controller_cmpc_2d/cost_functions']);

addpath([root, 'MPCQP/controller_predator']);

%% Local parameters
Runs = 1;
controller = 'cmpc';
experiment = 'bf';

%% Predator params
params.predator = 0;
params.pFactor = 1.40;
params.pred_radius = 6;
params.wp = 500;

%% Params
params.steps = 1000;
params.n = 2;
params.h = 6;
params.dt = 0.1;
params.ct = 0.2;
params.amax = 1.5; %1.5
params.vmax = 2;%2

params.minP = -5;
params.maxP = 5;
params.minV = 0;
params.maxV = 1;
params.Dmax = 100; % neighbor def
params.dmin = 2.0;
params.Ds = params.dmin; % safe dist

params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 0.05; %1e0; %nominal case
% params.rs = sensing_radius(params);
params.rs = 4;

params.ws = 39.06252;
params.wc = 1;
params.wo = 400;
params.wt = 10;

%% Load obstacle file.
rect_file = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/results_neurips/obstacles/3/rectangles.txt';
params.rects = readRects(rect_file);

run controller_cmpc_2d/common/create_mex.m;

%% Optimizer Settings
opt = optimoptions('fmincon');
opt.Display = 'off';
% opt.Algorithm = 'active-set';
opt.MaxIterations = 100;
% opt.MaxFunctionEvaluations = 6000;
% opt.StepTolerance = 1e-12;
% Stopping criteria if change in function value is lesser than this:
opt.FunctionTolerance = 1e-7;

%% Load initial states for multiple experiments
init_file_name = 'istate_15_1000';
load(['trajs_fossacs2020/src/' init_file_name]);

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
    %         [posi, veli] = gen_state_random(params);
    %         step_number = 1;
    %         load('trajs_fossacs2020/bf/cmpc/traj_cmpc_bf_1.mat');
    %         posi = [traj.x(step_number,:) ;traj.y(step_number,:)];
    %         veli = [traj.vx(step_number,:) ;traj.vy(step_number,:)];
    %     veli = zeros(size(veli));
    
    posi = [0.1, 3.04741; 0.0860029,  3.1];
    
    veli = [1, -0.996978; 1.0556, -1.09107];
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
    bd = zeros([params.n, params.steps, params.n]);
    mpc_cost = zeros(1, params.steps);
    f = zeros(1, params.steps);
    exit_flag_optimizer = zeros(1, params.steps);
    
    x(1,:) = posi(1,:);
    y(1,:) = posi(2,:);
    vx(1,:) = veli(1,:);
    vy(1,:) = veli(2,:);
    f(1) = fitness(posi, veli, params);
    
    
    %% Controller and Dynamics:
    pos = posi; vel = veli;
    prev_sol = zeros(2*params.n*params.h,1);
    for t = 1:params.steps
        if mod(t,5) == 0
            e = round(toc, 1);
            disp(['run: ' num2str(i) '/' num2str(Runs) ', step: ' num2str(t) '/' num2str(params.steps) ', Time:' num2str(e) 's, CPU:' num2str(cpu_time)]);
        end
        
        if mod(t - 1, controller_run) == 0
            if params.predator
                [acc, fval, e_flag, prev_sol, history] = controller_cmpc_2d(pos, vel, prev_sol, params, opt);
                
                [acc(:,end), dist_to_center] = controller_predator(pos, params);
                if dist_to_center < 10
                    acc(:,end) = [ax(t,end); ay(t,end)];
                end
            else
                start_cpu_time = cputime;
                [acc, fval, e_flag, prev_sol, history] = controller_cmpc_2d(pos, vel, params, opt);
                end_cpu_time = cputime;
                cpu_time = cpu_time + end_cpu_time - start_cpu_time;
                function_calls = function_calls + 1;
                %                  plot(history.fval, 'LineWidth', 1.5, 'Marker', '*', 'Color', 'r', 'MarkerEdgeColor', 'b');
                %                  pause(0.1);
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
        mpc_cost(t) = fval;
        exit_flag_optimizer(t) = e_flag;
        f(t+1) = fitness(pos, vel, params);
    end
    
    %% Store output.
    traj(i).x = x;
    traj(i).y = y;
    traj(i).vx = vx;
    traj(i).vy = vy;
    traj(i).ax = ax(1:params.steps,:);
    traj(i).ay = ay(1:params.steps,:);
    traj(i).mpc_cost = mpc_cost;
    traj(i).fitness = f;
    traj(i).exit_flag = exit_flag_optimizer;
    traj(i).params = params;
    
    displayTraj(x,y,vx,vy); title('Minimally Invasive', 'FontSize', 17);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, ['Images/trajComparison_' num2str(i) '.jpg']);
    pause(0.5);
end

%% Save output
out_traj_name = ['traj_' controller  '_' experiment '_' num2str(Runs)];
out_folder_name = ['trajs_fossacs2020/' experiment '/' controller  '/'];
save([out_folder_name out_traj_name], 'traj');
toc;


    










