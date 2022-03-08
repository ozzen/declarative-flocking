clc
clear all

%% dependencies
% root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
%
% addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-functions');
% addpath([root, 'MPCQP/common']);
%
% addpath([root, 'MPCQP/controller_cmpc_2d']);
% addpath([root, 'MPCQP/controller_cmpc_2d/common']);
% addpath([root, 'MPCQP/controller_cmpc_2d/cost_functions']);
%
% addpath([root, 'MPCQP/controller_predator']);

%% parpool options
% c = parcluster('local');
% c.NumWorkers = 32;
% parpool(32);

%% Relative paths

addpath('../m-functions');
addpath('Common');

addpath('controller_cmpc_2d');
addpath('controller_cmpc_2d/common');
addpath('controller_cmpc_2d/cost_functions');

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

params.ws = 39.06252;
params.wc = 1;
params.wo = 400;
params.wt = 10;

run controller_cmpc_2d/common/create_mex.m;

%% Generate states. set-up
samples_per_var = 12;
disp(['Original: ' num2str(samples_per_var^6)]);

%% pos
pos_per_dim = linspace(params.minP, params.maxP, samples_per_var);

pos = augment_matric(pos_per_dim');

%% vel
v_mag = linspace(params.vmax/samples_per_var,params.vmax, samples_per_var);
theta = linspace(0,360 * ((samples_per_var - 1) / samples_per_var), samples_per_var);

vx = v_mag' * cosd(theta);
vy = v_mag' * sind(theta);

vx = [0; vx(:)];
vy = [0; vy(:)];

vel = [vx, vy];

%% Name and save
save_to_dir = 'Training_data/';
file_name = ['samples_', num2str(params.n), '.txt'];

fid = fopen([save_to_dir file_name],'wt');
M = 'px1,px2,py1,py2,qx1,qx2,qy1,qy2,ax1,ax2,ay1,ay2\n';
fprintf(fid, M);

%% Nested loops approach
n_pos = size(pos,1);
n_vel = size(vel,1);
count = 0;
tic;
parfor p = 1:n_pos
    for v1 = 1:n_vel
        for v2 = 1:n_vel
            pos_k = [[0; 0] pos(p,:)'];
            vel_k = [vel(v1,:)' vel(v2,:)'];
            
            delta_p = pos(p,:);
            delta_v = vel(v2,:) - vel(v1,:);
            filter = barrier_function(delta_p, delta_v, params) >= 0;
            
            if filter
                [a, fval, e_flag, prev_sol, history] = controller_cmpc_2d(pos_k, vel_k, params, opt);
                M = [pos_k(:)' vel_k(:)' a(:)'];
                dlmwrite([save_to_dir file_name], M, '-append', 'delimiter', ',');
            end
            
            disp(['p:' num2str(p) ', v1:' num2str(v1) ', v2:' num2str(v2) ', states done: ' num2str((p-1)*n_vel*n_vel + (v1-1)*n_vel + v2)]);
            
        end
    end
end
toc;




