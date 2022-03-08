%% main_view_animation
% Display video for 3D flock for
%   1. 3D point model dynamics
%   2. Quadcopter model dynamics (6-DOF)
%
%
%%
clc
clear all
close all

%%
dest_path = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/Results/';
source_path = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';

%% dependencies
addpath([source_path 'm-scripts/']);
addpath([source_path 'm-functions/']);


%% Some on/off options
skip_1 = 5; %How many steps to skip in video between frames. 

Acceleration = 1;
Acc_scale = 1;
v_scale = 0.5;

isText = 0;
showQuad = 0;

drawTraj = 2; %0,1,2
trail_len = 400;


FullBox = 1;
BoxLength = 2;
BoxL = 1;

rotor_speed = 1;

%% load from single traj
file_name = '20_6_400_ [2020-06-03]'; 
folder_name = file_name;
fname = [dest_path, folder_name, '/', file_name, '.mat'];
load(fname, 'traj');

%% Load from the array of traj
% traj_name = 'cmpc_quad_bf.mat';
% folder_name = 'point_1';
% file_name = 'vid_2';
% load(['experiments_iccps2020/' traj_name], 'traj');
% traj = traj(11);
%%
[s,x,y,z,vx,vy,vz,accx,accy,accz,omegas,params] = read_data_mat(traj, showQuad);

%% find edge cluster
% params.num_leaders = 4;
% iter = 100;
% params.turn_angle = 90;
% [agent_ids] = edge_cluster(s(iter,:,:), params);
% [edge_agents, acc_turn, n] = edge_cluster(s(2*params.start_turn+1,:,:), params);

%% Draw trajecrory
f1 = figure;
for i = 1:params.n
    plot3( s(:,1,i), s(:,2,i), s(:,3,i) , 'LineWidth', 1.5, 'Color', [1,0,0]);
    hold on
    scatter3(s(1,1,i), s(1,2,i), s(1,3,i), 200, [0,1,0], '.');
end
axis equal

%% View Angle
% n = traj.axis_rot;
if params.turn
    [azi, ele] = vec2azi_ele(traj.axis_rot);
else
    azi = 35;
    ale = 35;
end

view(azi,ele)

%% Plot Normal
% vel = squeeze(s(params.start_turn+1,4:6,:));
% pos = squeeze(s(params.start_turn+1,1:3,:));
% 
% centroid = mean(pos,2);
% n = 5*n;
% plot3([centroid(1) n(1)+centroid(1)], [centroid(2) n(2)+centroid(2)], [centroid(3) n(3)+centroid(3)], 'k', 'LineWidth', 3)

%% Set Video options 
mkdir([dest_path folder_name])
videoName = [dest_path folder_name '/' file_name ];
myVideo = VideoWriter(videoName, 'MPEG-4');
fastForwardRate = 1;
myVideo.FrameRate = min(30,fastForwardRate * (1/params.dt));
myVideo.Quality = 50;
open(myVideo);
skip = skip_1 + 1;
%% Script to generate video
view_3D
