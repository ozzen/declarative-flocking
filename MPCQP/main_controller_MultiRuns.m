clc
clear
close all

% Controllers
addpath('controller_qp_constrained')
addpath('controller_mininvasive_centralized')

% From MPCQP
addpath('Common')
addpath('CostFunctions')
% From ..
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Li Wangs code/Multi_obj');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Plotting');

Runs = 1;
params.steps = 2000;
%% Params:
params.n = 2;
params.h = 4;
params.dt = 0.005;
params.ct = 0.005;
params.amax = 5;
params.vmax = 2.5;

params.minP = -8;
params.maxP = 8;
params.minV = 0;
params.maxV = 2;
params.Dmax = 100; % neigbor def
params.Ds = 1.0; % safe dist
params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 10; %1e0; %nominal case

options =  optimoptions('Display','off');

u = zeros(2*params.n*params.h + 1,1);
for i = 1:Runs
% %% Load init:
% %     [posi, veli] = readInit('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/experiments/exp1/init_conf_%d.txt', i);
% %     posi = [-8 8; -5 5];
    [posi, veli] = genInitConf(params);
%     
    veli = zeros(size(veli));
%     x = zeros([params.steps, params.n]);
%     y = zeros([params.steps, params.n]);
%     vx = zeros([params.steps, params.n]);
%     vy = zeros([params.steps, params.n]);
%     ax = zeros([params.steps, params.n]);
%     ay = zeros([params.steps, params.n]);
%     cost = zeros(1, params.steps);
% 
%     
% %% QP Controller and Dynamics:  
%     pos = posi; vel = veli;
%     for t = 1:params.steps
%         [ acc, fval, ~, ~, flag ] = controller_QP_constrained( params, pos, vel, u );
%         if flag, continue; end %Controller failed
%         [pos, vel] = Dynamics(pos, vel, acc, params);
%         x(t,:) = pos(1,:);
%         y(t,:) = pos(2,:);
%         vx(t,:) = vel(1,:);
%         vy(t,:) = vel(2,:);
%         ax(t,:) = acc(1,:);
%         ay(t,:) = acc(2,:);
%         cost(t) = fval;
%     end
    
%% Store output.
%     trajQPC(i).x = x;
%     trajQPC(i).y = y;
%     trajQPC(i).vx = vx;
%     trajQPC(i).vy = vy;
%     trajQPC(i).ax = ax;
%     trajQPC(i).ay = ay;
%     trajQPC(i).params = params;
%     
%     close all;
%     subplot(1,2,1);
%     displayTraj(x,y,vx,vy); title('MPC Constrained',  'FontSize', 17);

%% MI Controller and Dynamics:
    h_total = zeros(nchoosek(params.n, 2));
    pos = posi; vel = veli;
    for t = 1:params.steps
        [H, ~, ~, ~, Aeq, Beq, ~, lb, ub ] = GenMatQP( params, pos, vel, u);
        [u, fval, exitflag, output,lambda] = quadprog(H,[],[],[],Aeq,Beq,lb,ub,[],options);
         acc = u2acc( u, params );
        [acc, flag, ~, ~, h_allpairs, output2] = controller_min_invasive(acc, pos, vel, params);
        
        if flag, break; end %Controller failed

        h_total(:, t) = h_allpairs;
        [pos, vel] = Dynamics(pos, vel, acc, params);
        x(t,:) = pos(1,:);
        y(t,:) = pos(2,:);
        vx(t,:) = vel(1,:);
        vy(t,:) = vel(2,:);
        ax(t,:) = acc(1,:);
        ay(t,:) = acc(2,:);
        cost(t) = fval;
    end
    
%% Store output.
    trajMI(i).x = x;
    trajMI(i).y = y;
    trajMI(i).vx = vx;
    trajMI(i).vy = vy;
    trajMI(i).ax = ax;
    trajMI(i).ay = ay;
    trajMI(i).params = params;
    trajMI(i).h = h_total;
    
    subplot(1,2,2);
    displayTraj(x,y,vx,vy); 
    title('Minimally Invasive', 'FontSize', 17);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, ['Images/trajComparison_' num2str(i) '.jpg']);
    pause(0.5);
end

%Save results
% save('Results_QPConstrained.mat', 'trajQPC');
save('Results_MinInvasive.mat', 'trajMI');





