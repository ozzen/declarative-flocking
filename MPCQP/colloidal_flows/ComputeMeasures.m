clc
clear all
close all

% source = '/Users/Usama/Swarms New'; %DF
source = '/Users/Usama/reynolds/Reynolds Model C/'; %Reyn

load([source 'Results_Obst_0.mat']);
T0 = trajObst;
load([source 'Results_Obst_25.mat']);
T25 = trajObst;
load([source 'Results_Obst_50.mat']);
T50 = trajObst;
load([source 'Results_Obst_75.mat']);
T75 = trajObst;

params0 = T0.params;
params25 = T25.params;
params50 = T50.params;
params75 = T75.params;

%%
[ trajs_array_0 ] = AddMeasuresFn(T0);
[ trajs_array_25 ] = AddMeasuresFn(T25);
[ trajs_array_50 ] = AddMeasuresFn(T50);
[ trajs_array_75 ] = AddMeasuresFn(T75);


temp = trajs_array_50(1,1);
steps = size(temp{1}, 1);

%%
fnames = {'x' ,'y', 'vx', 'vy', 'ax', 'ay', 'params', 'minPW',...
    'maxPW', 'pnet', 'connectivity',  'v_converg', 'irreg', 'diameter'};

% plotsNeeded = [9, 11, 12, 13];
plotsNeeded = [14];
for fnum = plotsNeeded
    figure;
    h1 = printMeasure(fnames{fnum}, trajs_array_0, fnum, params0, 'r', '--');
    h2 = printMeasure(fnames{fnum}, trajs_array_25, fnum, params25, 'b', '--');
    h3 = printMeasure(fnames{fnum}, trajs_array_50, fnum, params50, 'k', '--');
    h4 = printMeasure(fnames{fnum}, trajs_array_75, fnum, params75, [0.2 0.4 0], ':');

    % plot([0 steps],[paramsMI.Ds paramsMI.Ds], 'g--');
    % axis([0 paramsMI.dt*steps 0 10]);

    K = gobjects(4,1);
    K(1) = h1;
    K(2) = h2;
    K(3) = h3;
    K(4) = h4;

    legend(K, {'Num of Obstacles = 0', 'Num of Obstacles = 25', 'Num of Obstacles = 50', 'Num of Obstacles = 75'},'FontSize',17);

    set(gcf, 'Position', 0.65* get(0, 'Screensize'))
    saveas(gcf, ['/Users/Usama/Swarms New/Obstacle Results/' fnames{fnum} '_Reyn_SD.jpg'])
end
    
    
    
    
    
   
    

