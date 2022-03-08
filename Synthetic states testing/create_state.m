clc
close all
clear

%%
root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
addpath([root 'm-functions/init_generators']);

addpath("/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_reynolds");
addpath("/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_baseline_barrier");

%% Params
params.n = 5;
params.h = 1;
params.dt = 0.05;
params.ct = 0.05;
params.amax = 5;
params.vmax = 2.5;
params.knn = 1;

params.minP = -5;
params.maxP = 5;
params.minV = -1;
params.maxV = 1;
params.Dmax = 100; % neighbor def

params.dmin = 1.0;
params.Ds = params.dmin; % safe dist

params.min = 2.0;
params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 0.05; %1e0; %nominal case
% params.rs = sensing_radius(params);
params.rs = 4;
params.rs_bc = 2;
% fsl_switch_boundary = switching_boundary(params, 8);
% rsl_switch_boundary = switching_boundary(params, 10);

fsl_switch_boundary = 1.50;
rsl_switch_boundary = 1.70;


weights.align = 0;%1
weights.cohes = 3;
weights.sep = 0;%2


%% optimizer settings
opt = optimoptions('fmincon');
opt.Display = 'off';
% opt.Algorithm = 'active-set';
opt.MaxIterations = 6000;
opt.MaxFunctionEvaluations = 6000;
opt.ConstraintTolerance = 1e-6;
opt.OptimalityTolerance = 1e-6;

% opt.StepTolerance = 1e-12;
% Stopping criteria if change in function value is lesser than this:
opt.FunctionTolerance = 1e-7;

%% set figure 

f = figure;
set(gcf, 'Position',  [400, 400, 1200, 1200])
axis equal
grid on
axis([-5, 5, -5, 5])
hold on

%% Generete State
% [pos, vel] = gen_state_manual(params, f);
% [pos, vel] = gen_state_radial(params, 3);
[pos, vel] = gen_state_random(params);

%% Display
scatter(pos(1,:), pos(2,:), 1000, 'r', '.');

pos_ = pos + vel;
plot([pos(1,:); pos_(1,:)], [pos(2,:); pos_(2,:)], 'b-');

%% Fixed state

% pos = [2.8565    1.7021; 4.2109    3.9897];
% vel = [-2.4974    2.0974; -0.1141    1.3605];

% displayInitState(pos(1,:), pos(2,:), vel(1,:), vel(2,:), 1, [1,0,0])

%% Save state

%% Run controller
[acc, fval] = controller_baseline_barrier(pos, vel, params, opt);

%% Barrier function
delta_p_ij = pos(:,1) - pos(:,2);
delta_v_ij = vel(:,1) - vel(:,2);

res = h_ij(delta_p_ij, delta_v_ij, [], [], params, 1);
title(num2str(res));

%% Display result
pos_ = pos(:,1) + acc';
plot([pos(1,1), pos_(1)], [pos(2,1), pos_(2)], 'g', 'LineWidth', 2);






