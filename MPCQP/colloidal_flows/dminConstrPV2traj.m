clc
clear all

%% Params:
params.n = 20;
params.h = 9;
params.dt = 0.3;
params.ct = 0.3;
params.amax = 3;
params.vmax = 9;

params.minP = -10;
params.maxP = 10;
params.minV = 0;
params.maxV = 2;
params.Dmax = 30; % neigbor def
params.Ds = 4.0; % safe dist
params.Dc = 0.8;%0.6; % connectivity dist
params.gamma = 10; %1e0; %nominal case

steps = 100;
x = zeros([steps, params.n]);
y = zeros([steps, params.n]);
vx = zeros([steps, params.n]);
vy = zeros([steps, params.n]);
ax = zeros([steps, params.n]);
ay = zeros([steps, params.n]);
%%
% NumObst = [50, 25, 0, 75];
NumObst = [0, 25, 50, 75];

for k = 1:4 
%     source = ['/Users/Usama/Swarms New/experiments/exp' num2str(k) '/gmpcPointObst75_%d_%s.txt'];
    source = ['/Users/Usama/reynolds/Reynolds Model C/experiments/exp' num2str(k) '/reynolds_%d_%s.txt'];
    
    runs = 20;

    for i = 1:runs
        [x, y, vx, vy] = read_data(i, params.n, source);

        trajObst(i).x = x;
        trajObst(i).y = y;
        trajObst(i).vx = vx;
        trajObst(i).vy = vy;
        trajObst(i).ax = ax;
        trajObst(i).ay = ay;
        trajObst(i).params = params;
    %     
    %     displayTraj(x,y,vx,vy); title('dMin Penalty',  'FontSize', 17);
    %     pause(3);
    %     close;
    end

    save(['Results_Obst_' num2str(NumObst(k)) '.mat'], 'trajObst');
end
