clc
clear
close all

%% Options and number of runs

numSim = 5;
addpath("/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Quadcopter_model_new");
dest_path = 'experiments_iccps2020/';

%% Params and codegen
params.n = 20;
params.h = 6;
params.dt = 0.05;
params.ct = 2 * params.dt;
params.vmax = 4;
params.amax = 2;

% Initialization options
params.ipos = [-3,3];
params.ivel = [1,2];


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
params.quad = 1;

%Prediction model [1:point, 2:quadcopter]
params.cmpc_prediction_model = 1;

% distributed
params.knn = 7;

%Create mex files
create_mex

%% optim options
opt = optimoptions('fmincon');
opt.Display = 'off';
% opt.Algorithm = 'active-set';
opt.MaxIterations = 3000;
opt.MaxFunctionEvaluations = 3000;
% opt.StepTolerance = 1e-12;
% Stopping criteria if change in function value is lesser than this:
opt.FunctionTolerance = 1e-2;
    
%% Init Gen.
% Generate
% init = gen_init(params, numSim);
% save([dest_path 'src_4000.mat'], 'init');
% Load
load([dest_path 'src_4000.mat']);

tic
%% Run Simulations
parfor i = 1:numSim
%     disp(['simulation number: ' num2str(i)]);
	[traj(i)] = RunSimulation(init{i}, params, i, opt);
end
toc
%% Save
save([dest_path 'cpmc_quad_bf_' num2str(numSim) '.mat'], 'traj');









