% clc
% clear
close all
addpath('/Users/Usama/Swarms New/Plotting');
addpath('Common');
addpath('CostFunctions');
addpath('/Users/Usama/Swarms New');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-functions');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/m-scripts');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Plotting')
addpath('controller_QP_constrained')
addpath('controller_mininvasive_centralized')
% addpath('/Users/Usama/Usama /Declarative Flocking/Li Wangs code/Multi_obj');
%%
params.n = 4;
params.h = 1;
params.dt = 0.05;
params.ct = 0.50;
% params.amax = 1.414;
% params.vmax = 11.314;
params.amax = 0.8;
params.vmax = 0.2;
steps = 1500;

wpControlVars;
params.minP = -10;
params.maxP = 10;
params.minV = 0;
params.maxV = 2;
params.Dmax = 30; % neigbor def
params.Ds = 0.2; % safe dist
params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 5e-1; %1e0; %nominal case
descrptn = 'wpMI';

options.Acc = 1;
options.makeVid = 1;
options.showVideo = true;

options.QPConstrained = 0;
options.MI = 1;
options.ConstDmin = 0;
options.wpControl = 0;
options.MIWayPoint = 0;

opt = optimoptions('quadprog');
opt.Display = 'off';


%% Gen Initial state.
[ posi, veli ] = genInitConf(params);
% [posi, veli] = readInit('/Users/Usama/Swarms New/experiments/exp1/init_conf_%d.txt', 1);

% posi = [0, 8; 5, 5];
% veli = [3, -3; 0.75, -0.75];

% posi = [4.0014    4.0023;    7.3975    2.6026];
% veli = [0.0226   -0.0223;   -0.2207    0.2209];
% sqrt(4*params.amax * (norm( posi(:,1) - posi(:,2) ) - params.Ds ) ); 
% posi = [-0.95, -0.45, 0.15, 0.65; 1.5, 1.5, 1.5, 1.5];

pos = posi;
vel = veli;
us = [];

cost = zeros(1, steps);
LieDerVec = zeros(1, steps);
alphaVec = zeros(1, steps);
u = zeros(2*params.n*params.h + 1, 1);
acc = zeros(size(vel));

x = zeros([steps, params.n]);
y = zeros([steps, params.n]);
vx = zeros([steps, params.n]);
vy = zeros([steps, params.n]);
ax = zeros([steps, params.n]);
ay = zeros([steps, params.n]);

x(1,:) = pos(1,:);
y(1,:) = pos(2,:);
vx(1,:) = vel(1,:);
vy(1,:) = vel(2,:);

%% Basic Plot
if options.showVideo
    i = 1;
    pos_ = pos + params.dt*vel;
    xD = [pos(1,:) ; pos_(1,:)];
    yD = [pos(2,:) ; pos_(2,:)];

    fg = figure;
    hold on
    axis equal
    grid on;
    CurAx = gca;
    CurAx.LineWidth = 0.5;
    CurAx.GridColor = [0.7 0.7 0.7];
    CurAx.GridAlpha = 1;
    hold on
    set(gcf, 'Position', get(0, 'Screensize'))

    c = repmat([1 0 0], params.n, 1);
    s = 400*ones(params.n,1);

    p = scatter(pos(1,:), pos(2,:), s, c ,'.');
    
    for ii = 1:params.n
        hold on
        l(ii) = quiver(pos(1,ii), pos(2,ii), vel(1,ii), vel(2,ii), 'b', 'linewidth', 1.5, 'showarrowhead', 'on','maxheadsize',0.5, 'autoscalefactor',0.3);
        ball(ii) = rectangle('position', [pos(:,ii)' - 0.5*[params.Ds params.Ds] [params.Ds params.Ds]], 'curvature', [1,1], 'LineWidth', 1.2, 'EdgeColor', 'b');
        if options.Acc
            ac(ii) = quiver(pos(1,ii), pos(2,ii), acc(1,ii), acc(2,ii), 'g', 'linewidth', 1.5, 'showarrowhead', 'on','maxheadsize',0.5, 'autoscalefactor',0.1);
        end
    end
%     l = plot(xD, yD, 'b', 'LineWidth', 1.3); %velocity plot
%     rect = rectangle('position', [pos(:,1)' - [params.Ds params.Ds] 2*params.Ds 2*params.Ds], 'Curvature', [1 1]);
    rect = 1;
    plot(wayPoints(:,1), wayPoints(:,2), '*', 'MarkerSize', 5);

%     if options.Acc
%         xA = [ pos(1,:); pos(1,:) + 10 * acc(1,:) ];
%         yA = [ pos(2,:); pos(2,:) + 10 * acc(2,:) ];
%         ac = plot(xA,yA,'g', 'LineWidth', 1); %acc plot
%     end
end

%% Directory setup
destPath = '/Users/Usama/Swarms New/MPCQP/Results/';
dirName = [num2str(params.n) '_' num2str(params.h) descrptn  ];
mkdir([destPath dirName]);
save([destPath dirName '/initState.mat'], 'posi', 'veli');

%% Set Video options 
videoName = [destPath dirName '/' dirName 'NoSol' ];
myVideo = VideoWriter(videoName, 'MPEG-4');
myVideo.FrameRate = 12;
myVideo.Quality = 100;
open(myVideo);

%% Set Window
if options.showVideo
    [axRange, limits] = setAxis( pos, params, [] );
    axis(axRange);
    textpos = (7/8) * mean(abs(limits));
    centre = sum(pos,2)/params.n;
    textt = text(centre(1), centre(2)+textpos, num2str(i));
end

strt = 2;
tic
for i = strt:steps
    t0 = toc;
    if options.showVideo
        currkey=get(gcf,'CurrentKey');
        if currkey == 'return'
            break;
        end
        if options.makeVid && mod(i , ceil(0.3/params.dt/myVideo.FrameRate)) == 0 
            F = getframe(fg);
            writeVideo(myVideo,F);
        end
    end
    %% Constrained MPC and Minimally Invasive
    if options.QPConstrained
        [ acc, fval, LD, alp, Flag ] = controller_QP_constrained( params, pos, vel, u );
        LieDerVec(i) = LD; 
        alphaVec(i) = alp;
        
    elseif options.MI
        [H, ~, ~, ~, Aeq, Beq, ~, lb, ub ] = GenMatQP( params, pos, vel, u);
        [u, fval, exitflag, output,lambda] = quadprog(H,[],[],[],Aeq,Beq,lb,ub,[],[],opt);
        %Back calculate acceleration.
        acc = u2acc( u, params );
        [ ~, velSat ] = Dynamics( pos, vel, acc, params );
        acc = (velSat - vel)./params.dt;
        %...........................
        [acc, flag, ME, LD, alp] = minInvasiveControl(acc, pos, vel, params);
        if (flag), break; end
        LieDerVec(i) = LD; 
        alphaVec(i) = alp;
 
    elseif options.ConstDmin
        [~, ~, ~, ~, ~, ~, ~, lb, ub ] = GenMatQP( params, pos, vel, u);
        [u, fval] = fmincon(@(u)costSum(u, pos, vel, params)...
            , zeros(2*params.n*params.h, 1), [],[],[],[],lb(1:end-1),ub(1:end-1), @(u)mycon(u, pos, vel, params));
%         [u, fval] = fmincon(@(u)costSum(u, pos, vel, params), zeros(2*params.n*params.h, 1), [],[],[],[],[],[]);

        %Back calculate acceleration.
        acc = u2acc( u, params );
        [ ~, velSat ] = Dynamics( pos, vel, acc, params );
        acc = (velSat - vel)./params.dt;
        %...........................
        
    elseif options.wpControl
        [~, ~, ~, ~, ~, ~, ~, lb, ub ] = GenMatQP( params, pos, vel, u);
        [u, fval] = fmincon(@(u)costSumWP( u, pos, vel, wayPoints, mde, target, params )...
            , zeros(2*params.n*params.h, 1), [],[],[],[],lb(1:end-1),ub(1:end-1), @(u)mycon(u, pos, vel, params));
        [ mde, Flg ] = updateMode( pos, wayPoints, mde, target, params );
%         if Flg
%             for ag = 1:params.n
%             plot(x(strt:2:i-1, ag), y(strt:2:i-1, ag), 'LineWidth', 1.4, 'Color', [0.5 0.5 1], 'LineStyle', '--');
%             end
%             strt = i;
%         end
        if mde > size(target,1), break; end
        %Back calculate acceleration.
        acc = u2acc( u, params );
        [ ~, velSat ] = Dynamics( pos, vel, acc, params );
        acc = (velSat - vel)./params.dt;
        
    elseif options.MIWayPoint
        [~, ~, ~, ~, ~, ~, ~, lb, ub ] = GenMatQP( params, pos, vel, u);
        [u, fval] = fmincon(@(u)costSumWP( u, pos, vel, wayPoints, mde, target, params )...
            , zeros(2*params.n*params.h, 1), [],[],[],[],lb(1:end-1),ub(1:end-1),  @(u)myConMI( u, params ));
        %Back calculate acceleration.
        acc = u2acc( u, params );
%         [ ~, velSat ] = Dynamics( pos, vel, acc, params );
%         acc = (velSat - vel)./params.dt;
        %...........................
        [acc, flag, ME, LD, alp] = minInvasiveControl(acc, pos, vel, params);
        [ mde, Flg ] = updateMode( pos, wayPoints, mde, target, params );
        if mde > size(target,1), break; end
        
        if (flag), break; end
        LieDerVec(i) = LD; 
        alphaVec(i) = alp;
    end
        
    %% Update 
    prevPos = pos;
    [pos, vel] = Dynamics(pos, vel, acc, params);
    x(i,:) = pos(1,:);
    y(i,:) = pos(2,:);
    vx(i,:) = vel(1,:);
    vy(i,:) = vel(2,:);
    ax(i,:) = acc(1,:);
    ay(i,:) = acc(2,:);
    cost(i) = fval;

    if options.showVideo
        setUpdates( ball, p, l, rect, ac, textt, pos, prevPos, vel, acc, limits, i, params, options);
    end
    t1 = toc;
    pause(params.dt/10);
end
steps = i;
% end

if options.makeVid
    F = getframe(fg);
    writeVideo(myVideo,F);
    close(myVideo);
end

%%
figure;
displayTraj( x(1:i-1,:), y(1:i-1,:) ,vx(1:i-1,:), vy(1:i-1,:));
saveas(gcf, [destPath dirName '/' 'traj' '.jpg']);


%%
figure;
i = 1;
rows = 4;
columns = 1;
%% Min Dist
subplot(rows, columns, i);
% try
%     openfig('fig2');
% catch
%     disp('fig not there.')
%     figure;
% end
K = PlotDistance(x, y, steps, params );
legend(K, {'Min pairwise distance', 'Max pairwise distance'},'FontSize',17);
i = i + 1;
%% BF
subplot(rows, columns, i);
[costBF] = PlotBf(x, y, vx, vy, params);
i = i + 1;

%% Fitness
subplot(rows, columns, i);
plotCurve(cost, 'Time', 'Cost', 1.5, params.n);
i = i + 1;

%% LD and alpha
subplot(rows, columns, i);
l1 = plotCurve(-real(LieDerVec), 'Time', 'LD, alpha', 1, params.n);
set(l1, 'Color', 'b', 'LineStyle', '-');
hold on 
l2 = plotCurve(-real(alphaVec), 'Time', 'LD, alpha', 1, params.n);
set(l2, 'Color', 'r', 'LineStyle', '-');
K1 = gobjects(2,1);
K1(1) = l1;
K1(2) = l2;
legend(K1, {'Lie-Derivative', '-alpha(B(x))'},'FontSize',17);
% axis([0 steps min(-real(LieDerVec)) max(-real(LieDerVec))]);
i = i + 1;

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, [destPath dirName '/' 'combined_plots1' '.jpg']);
savefig([destPath dirName '/' 'combined_plots1.fig']);


%% Other Graphs
figure;
%% Velocity plot
i = 1;
rows = 2;
cols = 2;
subplot(rows,cols,i);
plotVel( vx, vy );
i = i + 1;

%% Orientation plot
subplot(rows,cols,i);
PlotOrientation( vx(1:steps), vy(1:steps) );
i = i + 1;


%% Velocity fluctuation plot
subplot(rows,cols,i);
plotVelFluctuation( vx(1:steps), vy(1:steps));
i = i + 1;


%% Acceleration Plots and Fitness Plots
subplot(rows,cols,i);
plotAcc(vx, vy, params );
axis([0 steps 0 1.2 * params.amax])
i = i + 1;


%% Maximum Acc
% subplot(rows,cols,i);
% plotMaxAcc( vx, vy , params.n);
% i = i + 1;

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, [destPath dirName '/' 'combined_plots2' '.jpg']);
savefig([destPath dirName '/' 'combined_plots2.fig']);





