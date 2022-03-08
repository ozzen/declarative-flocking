clc
clear all
close all

addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Plotting');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-functions');
addpath('Common');
addpath('simplex_switching_logic');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/simplex_switching_logic_centralized');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter');

addpath('CostFunctions');
addpath('controller_qp_constrained')
addpath('controller_mininvasive_centralized')
addpath('controller_mininvasive_decentralized')
addpath('controller_reynolds')
addpath('controller_baseline_barrier')
addpath('controller_baseline_barrier_centralized')
addpath('controller_cmpc_2d')

addpath('plot_feasible');

addpath('Swarms New');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Li Wangs code/Multi_obj');

%% Params
Runs = 1;

params.steps = 20;
params.n = 15;
params.h = 6;
params.dt = 0.10;
params.ct = 0.20;
params.amax = 1.5;
params.vmax = 2;

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

% run controller_cmpc_2d/create_mex.m;
%%
fsl_switch_boundary = switching_boundary(params, 10);
rsl_switch_boundary = switching_boundary(params, 20);


weights.align = 0;%1
weights.cohes = 3;
weights.sep = 0;%2

u = zeros(2*params.n*params.h + 1,1);
indexes = 1:params.n;
acc = zeros(2, params.n);
controller_run = params.ct / params.dt;
zero_vec = zeros(params.n , 1) ;

opt = optimoptions('fmincon');
opt.Display = 'off';
opt.MaxIterations = 100;


for i = 1:Runs
    %% Load init:
    [posi, veli] = genInitConf(params);
%         step_number = 1;
%         load('results_minInvasive.mat');
%         posi = [traj.x(step_number,:) ;traj.y(step_number,:)];
%         veli = [traj.vx(step_number,:) ;traj.vy(step_number,:)];
%     
%     veli = zeros(size(veli));
    x = zeros([params.steps, params.n]);
    y = zeros([params.steps, params.n]);
    vx = zeros([params.steps, params.n]);
    vy = zeros([params.steps, params.n]);
    ax = zeros([params.steps+1, params.n]);
    ay = zeros([params.steps+1, params.n]);
    bd = zeros([params.n, params.steps, params.n]);
    policy = zeros([params.steps, 1]);
    cost = zeros(1, params.steps);
    mde = 1;
    
    %% MI Controller and Dynamics:
    pos = posi; vel = veli;
    
    for t = 1:params.steps
        if mod(t,5) == 0
            disp([ 'step ' num2str(t) ' out of ' num2str(params.steps)]);
        end
        
        if mod(t - 1, controller_run) == 0
            % Switching logic
%             if mde == 1
%                 if FSL_centralized(pos, vel, params, fsl_switch_boundary)
%                     mde = 2;
%                 end
%             elseif mde == 2
%                 if RSL_centralized(pos, vel, params, rsl_switch_boundary)
%                     mde = 1;
%                 end
%             end
            % control
            if mde == 1
%                 % Advanced controller (returns acc)
%                 for ai = 1:params.n
%                     %Neighbors
%                     distance_to_ai = sqrt(sum((pos - pos(:,ai)).^2, 1));
%                     neighbors = find(distance_to_ai < params.rs & distance_to_ai > 0 );
% 
%                     %Reynolds Control
%                     posN = [pos(:,ai), pos(:, neighbors)];
%                     velN = [vel(:,ai), vel(:, neighbors)];
%                     
%                     % reynolds's model
%                     acc(:,ai) = controller_reynolds(posN, velN, params, weights);
%                 end
                % cmpc 
                acc = controller_cmpc_2d(pos, vel, params, opt);
            else 
                % Baseline controller (returns acc)
                [acc, fval] = controller_baseline_barrier_centralized(pos, vel, params);
                cost(t) = fval;
            end
        else
            acc = [ax(t, :); ay(t, :)]; 
        end
        
        %h_total(:, t) = h_allpairs;
        [pos, vel] = Dynamics(pos, vel, acc, params);
        x(t,:) = pos(1,:);
        y(t,:) = pos(2,:);
        vx(t,:) = vel(1,:);
        vy(t,:) = vel(2,:);
        ax(t+1,:) = acc(1,:);
        ay(t+1,:) = acc(2,:);
        policy(t,:) = mde;
        %         linear_A{t} = A;
        %         linear_b{t} = b;
    end
    
    %% Store output.
    traj(i).x = x;
    traj(i).y = y;
    traj(i).vx = vx;
    traj(i).vy = vy;
    traj(i).ax = ax(1:params.steps,:);
    traj(i).ay = ay(1:params.steps,:);
    traj(i).params = params;
    traj(i).policy = policy;
    traj(i).bd = bd;
    %     traj(i).linear_A = linear_A;
    %     traj(i).linear_b = linear_b;
    
%     subplot(1,2,1);
    displayTraj(x,y,vx,vy); title('Minimally Invasive', 'FontSize', 17);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, ['Images/trajComparison_' num2str(i) '.jpg']);
    pause(0.5);
end
save('traj_cmpc.mat', 'traj');








