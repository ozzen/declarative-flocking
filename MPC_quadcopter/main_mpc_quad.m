%% main_mpc_quad - script used for 3D point-model and Quadrotor experiments.
clc
clear
close all

%% Dependencies
root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
addpath([root 'Quadcopter_model_new']);

addpath([root 'MPC_quadcopter/controller_cmpc_3d']);
addpath([root 'MPC_quadcopter/controller_cmpc_3d/common']);
addpath([root 'MPC_quadcopter/controller_cmpc_3d/cost_functions']);

addpath([root 'MPC_quadcopter/controller_dmpc_3d']);
addpath([root 'MPC_quadcopter/controller_dmpc_3d/common']);
addpath([root 'MPC_quadcopter/controller_dmpc_3d/cost_functions']);

%% Params and codegen
params.n = 20;
params.h = 6;
params.dt = 0.05;
params.ct = 2 * params.dt;
params.vmax = 4;
params.amax = 2;

% Initialization options
params.ipos = [-3,3];
params.ivel = [0,1];


params.delta_angle = deg2rad(30); % not used

params.t_end = params.dt * 400;
params.steps = params.t_end / params.dt;
params.m =0.650;
params.L = 0.23;
params.dmin = 4 * params.L;
params.g = 9.81;
params.j_r = 6e-5;
params.Ixx = 7.5e-3;
params.Iyy = 7.5e-3;
params.Izz = 1.3e-2;

%bounds
params.max_rotor_speed = inf; %rad/sec
params.max_angular_speed = 10; %rad/sec
params.max_angular_acc = 10; %rad/sec2

%angular speed to thrust
params.k_f = 3.13e-5; %b
params.k_m = 7.5e-7; %d

%wfitness
params.wc = 5;
params.ws = 30;

% plant model -> 1:Quadrotor, 0:point-model
params.quad = 0;

%Prediction model [1:point, 2:quadcopter]
params.cmpc_prediction_model = 1;

% distributed
params.knn = 7;

%Flock Turning
params.turn = 1; %0:off, 1:on
params.eta = 400;
params.w_m = zeros(1,params.n);

params.start_turn = 200;
params.num_leaders = 5;
params.t_fix = 70;

params.turn_angle = -90;
params.turn_acc = 0.25 * params.amax;


%Create mex files
run controller_dmpc_3d/common/create_mex.m;

%% optim options
opt = optimoptions('fmincon');
opt.Display = 'off';
% opt.Algorithm = 'active-set';
opt.MaxIterations = 3000;
opt.MaxFunctionEvaluations = 3000;
% opt.StepTolerance = 1e-12;
% Stopping criteria if change in function value is lesser than this:
opt.FunctionTolerance = 1e-2;
    
%%
destPath = [root 'MPC_quadcopter/Results/'];
naming_convention = [params.n, params.h, params.steps];
date_string = datestr(datetime,' [yyyy-mm-dd]');
fname = strcat(strjoin(string(naming_convention),'_'), '_', date_string);
mkdir(strcat(destPath, fname));

%% init config
pos = params.ipos(1) + (params.ipos(2)-params.ipos(1)) * rand(1,3, params.n);
vel = params.ivel(1) + (params.ivel(2)-params.ivel(1)) * rand(1,3, params.n);
% vel = zeros(1,3, params.n);
omegas = zeros(1,3,params.n);
% for i = 1:params.n
%     v = squeeze(vel(:,:,i))';
%     angles = vec2euler(v);
%     omegas(1,:,i) = [angles(3), angles(2), angles(1)];
% end
rest = zeros(1,3,params.n);
s_init = cat(2, pos, vel, omegas, rest);

%% Run Simulation

[traj] = RunSimulation(s_init, params, 1, opt);
s = traj.s;
accs = traj.accs;

%% Plot trajectories

f1 = figure;
for i = 1:params.n
    plot3( s(:,1,i), s(:,2,i), s(:,3,i) , 'LineWidth', 1.5, 'Color', [1,0,0]);
    hold on
    scatter3(s(1,1,i), s(1,2,i), s(1,3,i), 200, [0,1,0], '.');
end
axis equal

saveas(gcf, strcat(destPath, fname, '/', fname,'_path.jpg'));
savefig(strcat(destPath, fname, '/', fname, '_path.fig'));

%% Plotting
f2 = figure;
set(f2, 'Position', [600,100,1200,1200])
num_plots = 13;
t = 0:params.dt:params.t_end-params.dt;
for i = 1:params.n
    subplot(num_plots,1,1);
    plot_wrapper(t, s(1:params.steps,1,i), 1.5, '-', [1,0,0], '', 'x');
    hold on;
    
    subplot(num_plots,1,2);
    plot_wrapper(t, s(1:params.steps,2,i), 1.5, '-', [1,0,0], '', 'y');
    hold on
    
    subplot(num_plots,1,3);
    plot_wrapper(t, s(1:params.steps,3,i), 1.5, '-', [1,0,0], '', 'z');
    hold on
    
    subplot(num_plots,1,4);
    plot_wrapper(t, s(1:params.steps,4,i), 1.5, '-', [0,1,0], '', '$v_x$');
    hold on
    
    subplot(num_plots,1,5);
    plot_wrapper(t, s(1:params.steps,5,i), 1.5, '-', [0,1,0], '', '$v_y$');
    hold on
    
    subplot(num_plots,1,6);
    plot_wrapper(t, s(1:params.steps,6,i), 1.5, '-', [0,1,0], '', '$v_z$');
    hold on
    
    subplot(num_plots,1,7);
    plot_wrapper(t, s(1:params.steps,7,i), 1.5, '-', [1,0,0], '', '$\phi$/x');
    hold on
    
    subplot(num_plots,1,8);
    plot_wrapper(t, s(1:params.steps,8,i), 1.5, '-', [1,0,0], '', '$\theta$/y');
    hold on
    
    subplot(num_plots,1,9);
    plot_wrapper(t, s(1:params.steps,9,i), 1.5, '-', [1,0,0], '', '$\psi$/z');
    hold on
    
    subplot(num_plots,1,10);
    plot_wrapper(t, accs(1:params.steps,1,i), 1.5, '-', [0,0,1], '', '$a_x$');
    hold on
    
    subplot(num_plots,1,11);
    plot_wrapper(t, accs(1:params.steps,2,i), 1.5, '-', [0,0,1], '', '$a_y$');
    hold on
    
    subplot(num_plots,1,12);
    plot_wrapper(t, accs(1:params.steps,3,i), 1.5, '-', [0,0,1], '', '$a_z$');
    hold on
    subplot(num_plots,1,13);
    speed = sqrt(sum(s(1:params.steps,4:6,i).^2, 2))';
    plot_wrapper(t, speed , 1.5, '-', [1,0,1], '', '$speed$');
    hold on
end

saveas(gcf, strcat(destPath, fname, '/',fname, '_s.jpg'));
savefig(strcat(destPath, fname, '/',fname, '_s.fig'));

%% Plot fitness
f3 = figure;
subplot(2,1,1)
plot_wrapper(t, traj.mpc_cost, 1.5, '-', [1,0,0], '', 'fitness');

subplot(2,1,2)
plot_wrapper(t, traj.exit_flag, 1.5, '-', [1,0,0], '', 'e_flag');

saveas(gcf, strcat(destPath, fname, '/',fname, '_fit.jpg'));
savefig(strcat(destPath, fname, '/',fname, '_fit.fig'));

%% save trajectory params.steps x params.n
save(strcat(destPath, fname, '/', fname, '.mat'), 'traj');









