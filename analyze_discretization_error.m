clc
clear all
source_path = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
addpath([source_path 'Plotting/']);
addpath([source_path 'm-scripts/']);
addpath([source_path 'm-functions/']);
addpath([source_path 'performance_measures/'])
%%
src_folder = 'exp15';
prefix = ['experiments_ml/' src_folder '/'] ; %Source path
confFile = [prefix 'result_files/'  'gmpc.conf'];
[params ] = readConf(confFile);

pos_i = [0,0];
vel_i = [params.vmax/4,0];
acc = [0,params.amax];

params.steps = 20;

pos(:,1) = pos_i;
vel(:,1) = vel_i;

T = 1.02; 
f1 = figure;
f2 = figure;
for dt = 0.01:0.02:0.20
    params.steps = uint8(T/dt);
    pos = zeros(2,params.steps);
    vel = zeros(2,params.steps);
    pos(:,1) = pos_i;
    vel(:,1) = vel_i;
    for i = 2:uint8(T/dt)
        vel(:,i) = vel(:,i-1) + dt * acc';
        if norm(vel(:,i)) > params.vmax
            vel(:,i) = vel(:,i) * (params.vmax / norm(vel(:,i)) );
        end
        pos(:,i) = pos(:,i-1) + dt * vel(:,i);
    end
    figure(f2)
    plot(dt*double((1:uint8(T/dt))), sqrt(sum(vel.^2,1)), '.-');
    axis equal 
    hold on

    figure(f1);
    plot(pos(1,1:end-1), pos(2,1:end-1),'.-');
    axis equal
    hold on
end
