clc
clear all
close all

addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Plotting');
addpath('Common');
addpath('CostFunctions');
addpath('Swarms New');
addpath('/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Li Wangs code/Multi_obj');
%%
params.n = 4;
params.dt = 0.05;
params.steps = 150;
params.amax = 11;
params.vmax = 1.4;

p_circ = [-0.95,1.1; -0.45,1.1; 0.15,1.1; 0.65,1.1; 
          -0.95,0.4; -0.45,0.4; 0.15,0.4; 0.65,0.4;
          -0.95,-0.3; -0.45,-0.3; 0.15,-0.3; 0.65,-0.3;
          -0.95,-1; -0.45,-1; 0.15,-1; 0.65,-1];%goal positions
      
gd = 1+[0,1,2,3;7,6,5,4;8,10,9,11;15,13,12,14]; %goal index

gain = 1:1/params.n:2;
posi = zeros(params.n, 2);
veli = zeros(params.n, 2);
uhat = zeros(params.n, 2);
fdone = zeros(1,params.n);
fk = 1 ;
%%
pos = posi;
vel = veli;
velSat = veli;
for i = 1:params.steps
    if i==1, u = uhat; end %initialize current velocity
    for j = 1:params.n
        [uhat(j,:), fdone(j)] = GoTo_2nd_order_New(pos(j,:), vel(j,:),p_circ(gd(fk,j),:), gain(j),params.amax);
        velSat(j,:) = vel(j,:) + uhat(j,:)*params.dt;
        if norm(velSat(j,:)) > params.vmax
            velSat(j,:) = velSat(j,:)/norm(velSat(j,:))*params.vmax; %constraint max velocity
        end
        if sum(fdone)== params.n && fk<4
            fk = fk+1 ;
        end
    end
    %saturatioin on velocity
    uhat = (velSat - vel)/params.dt; % back calculate saturated acceleration
end