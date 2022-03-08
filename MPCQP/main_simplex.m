clc
clear
close all

root = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
addpath([root 'Plotting']);
addpath([root 'm-functions']);
addpath([root 'm-functions/init_generators']);

addpath([root 'performance_measures']);
addpath('Common');
addpath('simplex_switching_logic');
addpath([root 'MPC_quadcopter']);

addpath('CostFunctions');
addpath('controller_qp_constrained')
addpath('controller_mininvasive_centralized')
addpath('controller_mininvasive_decentralized')
addpath('controller_reynolds')
addpath('controller_baseline_barrier')  
addpath('controller_dnn_decentralized')
addpath('plot_feasible_region');

addpath('Swarms New');
addpath([root 'Li Wangs code/Multi_obj']);

%%
params.predator = 0;
params.pFactor = 1.40;
params.pred_radius = 6;
params.wp = 500;

%% Params for 2D neural flocking
Runs = 1;

params.steps = 500;
params.n = 15;
params.h = 1;
params.dt = 0.1; %0.05
params.ct = 0.1; % 0.05s
params.amax = 5;
params.vmax = 2.5;
params.knn = 1;

params.minP = -5;
params.maxP = 5 ;
params.minV = -1;
params.maxV = 1;
params.Dmax = 100; % neighbor def

params.dmin = 2.0;
params.Ds = params.dmin; % safe dist

params.min = 2.0;
params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 1.0; %1e0; %nominal case
% params.rs = sensing_radius(params);
params.rs = 4;
params.rs_bc = 4;
% params.fsl_switch_boundary = switching_boundary(params, 8);
% params.rsl_switch_boundary = switching_boundary(params, 10);

params.fsl_switch_boundary = 1.2; %TODO add this to params for record keeping %1.5
params.rsl_switch_boundary = 1.6; %1.9 


weights.align = 0;%1
weights.cohes = 3;
weights.sep = 0;%2

%%
u = zeros(2*params.n*params.h + 1,1);
indexes = 1:params.n;
acc = zeros(2, params.n);
controller_run = params.ct / params.dt;

%% Load NN
% modelfile = [root 'MPCQP/controller_dnn_decentralized/model_2agent.h5'];
% net = importKerasNetwork(modelfile);

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


%%
tic;
for i = 1:Runs
    %% Load init:
%     [posi, veli] = gen_state_random(params);
%     [posi, veli] = gen_state_radial(params, 3);
    [posi, veli] = gen_state_random_conditional(params, params.fsl_switch_boundary);
    
%     step_number = 1;
%     sim_num = 17;
%     load('traj_sim.mat');
%     posi = [traj(sim_num).x(step_number,:) ;traj(sim_num).y(step_number,:)];
%     veli = [traj(sim_num).vx(step_number,:) ;traj(sim_num).vy(step_number,:)];
    
    x = zeros([params.steps, params.n]);
    y = zeros([params.steps, params.n]);
    vx = zeros([params.steps, params.n]);
    vy = zeros([params.steps, params.n]);
    ax = zeros([params.steps+1, params.n]);
    ay = zeros([params.steps+1, params.n]);
    
    x(1,:) = posi(1,:);
    y(1,:) = posi(2,:);
    vx(1,:) = veli(1,:);
    vy(1,:) = veli(2,:);
    
    bd = zeros([params.n, params.steps, params.n]);
    h_ij = zeros([params.n, params.steps, params.n]);
    policy = zeros([params.steps, params.n]);
    cost = zeros(1, params.steps);
    mde = ones(1,params.n);
    mpc_cost = zeros(1, params.steps);
    %     f = zeros(1, params.steps);
    exit_flag_optimizer = zeros(1, params.steps);
    
    %% MI Controller and Dynamics:
    h_total = zeros(nchoosek(params.n, 2));
    pos = posi; vel = veli;
%     fig_convex = figure;
    fval = 0;
    e_flag = 0;
    %     hold on
    %     axis equal
    for t = 1:params.steps
        if mod(t,5) == 0
            e = round(toc, 1);
            disp(['run: ' num2str(i) '/' num2str(Runs) ', step: ' num2str(t) '/' num2str(params.steps) ', Time:' num2str(e) 's']);
        end
        if mod(t - 1, controller_run) == 0
            for ai = 1:params.n
                
                if mde(ai) == 1
                    [posN, velN, neighbors] = get_neighbors_states(ai, pos, vel, params.rs);
                    if FSL(posN, velN, params)
                        mde(ai) = 2;
                        [posN, velN, neighbors] = get_neighbors_states(ai, pos, vel, params.rs_bc);
                    end
                elseif mde(ai) == 2
                    [posN, velN, neighbors] = get_neighbors_states(ai, pos, vel, params.rs_bc);
                    if RSL(posN, velN, params)
                        mde(ai) = 1;
                        [posN, velN, neighbors] = get_neighbors_states(ai, pos, vel, params.rs);
                    end
                end
                
                if mde(ai) == 1
                    acci = controller_reynolds(posN, velN, params, weights);
                    % Deep Neural controller
                    %                     acci = controller_dnn_2d(pos ,vel, ai, params, net);
                elseif mde(ai) == 2
                    %[acci_new, flag, ME, h_allpairs, output, A, b] = controller_mininvasive_dec(acci,posN,velN,params);
                    [acci, fval, e_flag, history] = controller_baseline_barrier(posN, velN, params, opt);
%                     plot(history.fval, 'LineWidth', 1.5, 'Marker', '*', 'Color', 'r', 'MarkerEdgeColor', 'b');
%                     pause(0.1);
                    %                      acci = acci(end,:)';
                    %                     if ai == 1
                    %                         bd(:, t) = break_down;
                    %                     end
                end
                
                acc(:, ai) = acci;
                %% Plot convex region for one agent
                
                %                 if ai == 1
                %                     [l, f] = plot_convex(fig_convex, ai, pos, vel, params);
                %                     axis equal
                %                     axis([-5,5,-5,5]);
                %
                %                     pause(0.5);
                %                     delete(l);
                %                     delete(f);
                %                 end
                
%% store Lie Derivative & barrier funtion for all pairs.

                zero_vec_ld = zeros(params.n , 1);
                zero_vec_bf = zeros(params.n , 1) ;

                [~, lie_der, bar_fun] = cost_lie(acci, posN, velN, params);
                
                zero_vec_ld(neighbors) = lie_der;
                bd(:, t, ai) =  zero_vec_ld;
                
                zero_vec_bf(neighbors) = bar_fun;
                h_ij(:,t,ai) =  zero_vec_bf;

             end
        else
            acc = [ax(t, :); ay(t, :)];
        end
        
        %h_total(:, t) = h_allpairs;
        [pos, vel] = Dynamics(pos, vel, acc, params);
        x(t+1,:) = pos(1,:);
        y(t+1,:) = pos(2,:);
        vx(t+1,:) = vel(1,:);
        vy(t+1,:) = vel(2,:);
        ax(t+1,:) = acc(1,:);
        ay(t+1,:) = acc(2,:);
        policy(t,:) = mde;
        mpc_cost(t) = fval;
        exit_flag_optimizer(t) = e_flag;
        %         f(t) = fitness(pos, vel, params);
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
    traj(i).bd = bd; % b-derivative aka lie-derivative
    traj(i).h_ij = h_ij; % barrier function
    
    traj(i).mpc_cost = mpc_cost;
%     traj(i).fitness = f;
    traj(i).exit_flag = exit_flag_optimizer;
    %     traj(i).linear_A = linear_A;
    %     traj(i).linear_b = linear_b;
    
    %subplot(1,2,1);
    displayTraj(x,y,vx,vy); title('Minimally Invasive', 'FontSize', 17);
%     set(gcf, 'Position', get(0, 'Screensize'));
    %     pause(0.5);
end

%% Save output

date_string = datestr(datetime,' [yyyy-mm-dd]');
out_traj_name = ['traj_sim_' num2str(Runs)  date_string];

dir_path = [root 'MPCQP/traj/simplex/'];
% dest_path = [dir_path out_traj_name];

[dest_path, out_traj_name] = create_dir(dir_path, out_traj_name);
mkdir([dest_path '/Results']);

fprintf(['\nOUTPUT:\n' dest_path '\n']);

save([dest_path '/' out_traj_name '.mat'], 'traj');

% saveas(gcf, 'Images/traj_simplex.jpg');
toc;







