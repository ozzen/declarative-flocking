%% initial setup
directory = [destPath dirName '/'] ;
fitnessFile = [prefix '/gmpcOBST_%d_fitness.txt'];
fitnessFile = sprintf(fitnessFile, simNum);
LDFile = [prefix '/gmpcOBST_%d_LD.txt'];
LDFile = sprintf(LDFile, simNum);

%% Velocity plot
plotVel( vx, vy );
saveas(gcf, [directory 'velocity' '.jpg'])

%% Orientation plot
PlotOrientation( vx, vy )
saveas(gcf, [directory 'Orientation' '.jpg'])

%% Velocity fluctuation plot
plotVelFluctuation( vx, vy)
saveas(gcf, [directory 'Fluctuation in vel' '.jpg'])

%% Acceleration Plots and Fitness Plots
plotAcc(vx, vy )
saveas(gcf, [directory 'acc' '.jpg'])

%% Maximum Acc
plotMaxAcc( vx, vy , Num)
saveas(gcf, [directory 'max_acc' '.jpg'])

%% Fitness
plotFileData(fitnessFile, 'Fitness Value' ,Num);
saveas(gcf, [directory 'fit' '.jpg'])

%% LD
plotFileData(LDFile, 'Lie-Derivative' ,Num);
saveas(gcf, [directory 'LD' '.jpg'])

%% Display Traj
% subplot(1,2,1)
% displayInitState( x, y, vx, vy, 1)
% axis equal
% axis([-15 15 -15 15])

% subplot(1,2,2)
displayTraj( x, y ,vx, vy);
saveas(gcf, [directory 'Traj' '.jpg'])
%% Min Dist
figure;
PlotDistance;
set(gcf, 'Position', get(0, 'Screensize')./2)
legend(H, {'Min pairwise distance', 'Max pairwise distance'},'FontSize',17);
saveas(gcf, [directory 'minDist' '.jpg'])

close all;