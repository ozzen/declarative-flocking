clc
clear all

%% dependencies
root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';

addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-functions');
addpath([root, 'MPCQP/common']);

addpath([root, 'MPCQP/controller_cmpc_2d']);
addpath([root, 'MPCQP/controller_cmpc_2d/common']);
addpath([root, 'MPCQP/controller_cmpc_2d/cost_functions']);

addpath([root, 'MPCQP/controller_predator']);

%% Optimizer Settings
opt = optimoptions('fmincon');
opt.Display = 'off';
% opt.Algorithm = 'active-set';
opt.MaxIterations = 100;
% opt.MaxFunctionEvaluations = 6000;
% opt.StepTolerance = 1e-12;
% Stopping criteria if change in function value is lesser than this:
opt.FunctionTolerance = 1e-7;
%%
params.predator = 0;
params.pFactor = 1.40;
params.pred_radius = 6;
params.wp = 500;

%% Params 
params.steps = 36; % Number of samples
params.n = 2;
params.h = 6;
params.dt = 0.10;
params.ct = 0.30;
params.amax = 1.5; %1.5
params.vmax = 2;%2

params.minP = -5;
params.maxP = 5;
params.minV = -1;
params.maxV = 1;
params.Dmax = 100; % neighbor def
params.dmin = 2.0;
params.Ds = params.dmin; % safe dist

params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 0.05; %1e0; %nominal case
% params.rs = sensing_radius(params);
params.rs = 4;

params.ws = 440;
params.wc = 1;
params.wo = 400;
params.wt = 10;

run controller_cmpc_2d/common/create_mex.m;

%% Generate states.
samples_per_var = 2;
disp(['Original: ' num2str(samples_per_var^6)]);

v_mag = linspace(params.vmax/samples_per_var,params.vmax, samples_per_var);
theta = linspace(0,360, samples_per_var);

pos_per_dim = linspace(params.minP, params.maxP, samples_per_var);

vx = v_mag' * cosd(theta);
vy = v_mag' * sind(theta);

vx = vx(:);
vy = vy(:);

vel = augment_matric([vx, vy]);

pos = augment_matric(pos_per_dim');

% filter, based on distance
pos = pos(sum(pos.^2,2) >= params.dmin^2,:);

% data = [repmat(pos, samples_per_var^2, 1), vel];
data = augment_matric_pair(pos, vel);
data = [zeros(size(data, 1), 2), data];
num_data_points = size(data, 1);
disp(['After pos filter: ' num2str(num_data_points)]);

%% filter, based on barrier functions
delta_p = data(:,3:4);
delta_v = data(:,7:8) - data(:,5:6);
fltr = false(num_data_points, 1);
for i = 1:num_data_points
    fltr(i) = barrier_function(delta_p(i,:), delta_v(i,:), params) >= 0;
end

data = data(fltr, :);
num_data_points = sum(fltr);
disp(['After barrier filter: ' num2str(num_data_points)]);

computation_time = num_data_points * 0.31 / 3600;
disp(['computation_hours: ', num2str(computation_time)]);

%% Store in file.
save_to_dir = 'Training_data/';
out_file_name = ['halfSamples_', num2str(params.n), '_', num2str(num_data_points) '.txt'];
dlmwrite([save_to_dir out_file_name], data, ' ');

%% Run controller on states & store results.
tic;
q = zeros(num_data_points, params.n, 2);
p = zeros(num_data_points, params.n, 2);
a = zeros(num_data_points, params.n, 2);
% tic;
elapsed_time = 0;
for k = 1:num_data_points
    pos_k = [data(k,1:2); data(k,3:4)];
    vel_k = [data(k,5:6); data(k,7:8)];
    
    tic;
    [a(k,:,:), fval, e_flag, prev_sol, history] = controller_cmpc_2d(pos_k, vel_k, params, opt);
    one_run_time = toc;
    elapsed_time = elapsed_time + one_run_time;
%    	plot(history.fval, 'LineWidth', 1.5, 'Marker', '*', 'Color', 'r', 'MarkerEdgeColor', 'b');
%     pause(0.1);
    

    q(k,:,:) = pos_k;
    p(k,:,:) = vel_k;
    
    if mod(k-1,5) == 0
        e = round(elapsed_time, 1);
        disp(['states done: ' num2str(k-1) ', Time:' num2str(e) 's']);
    end
end
toc;
traj.q = q;
traj.p = p;
traj.a = a;
traj.params = params;

%% Name and save
save_to_dir = 'Training_data/';
traj_name = ['traj_', num2str(params.n), '_', num2str(num_data_points) '.mat'];
save([save_to_dir traj_name], 'traj');



