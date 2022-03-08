clc
clear all
close all
%% dirs
dir = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/figures/';
source_trajs = 'experiments_iccps2020/';


%% Loading trajs

traj_names = {'dnc_quad_bf_noscale.mat', 'dnc_point_bf_noscale.mat', 'cmpc_quad_bf.mat', 'cmpc_point_bf.mat'};

temp = load([source_trajs traj_names{1}]);
dnc_quad_traj = temp.traj;
% dnc_quad = squeeze(struct2cell(dnc_quad_traj));

temp = load([source_trajs traj_names{2}]);
dnc_point_traj = temp.traj;
% dnc_point = squeeze(struct2cell(dnc_point_traj));

temp = load([source_trajs traj_names{3}]);
cmpc_quad_traj = temp.traj;
% cmpc_quad = squeeze(struct2cell(cmpc_quad_traj));

temp = load([source_trajs traj_names{4}]);
cmpc_point_traj = temp.traj;
% cmpc_point = squeeze(struct2cell(cmpc_point_traj));

%% 1 - diameter 22,14,21,13
figure;
figName = 'diameter_diff';

t = linspace(0,20,400);
delta_dnc_diameter = mean_of_field(dnc_quad_traj, 22) - mean_of_field(dnc_point_traj, 14);
delta_cmpc_diameter = mean_of_field(cmpc_quad_traj, 21) - mean_of_field(cmpc_point_traj, 13);

% subplot(2,1,1)
plot(t,delta_dnc_diameter, 'r', 'LineWidth', 1.5);
hold on
plot(t,delta_cmpc_diameter, 'b', 'LineWidth', 1.5);
plot([0,20], [0,0], 'k');
ylabel('\Delta D', 'FontSize', 17);
xlabel('Time/s');
% legend('dnc-quad - dnc-point', 'cmpc-quad - cmpc-point', 'FontSize', 17)
legend('DNC-DNN', 'CMPC')
savefig(gcf, [dir figName '.fig']);


%% 2 - veloclity convergence 25,17,24,16
figure;
figName = 'vc_diff';

delta_dnc_vc = mean_of_field(dnc_quad_traj, 25) - mean_of_field(dnc_point_traj, 17);
delta_cmpc_vc = mean_of_field(cmpc_quad_traj, 24) - mean_of_field(cmpc_point_traj, 16);

% subplot(2,1,2)
plot(t,delta_dnc_vc, 'r', 'LineWidth', 1.5);
hold on
plot(t,delta_cmpc_vc, 'b', 'LineWidth', 1.5);
plot([0,20], [0,0], 'k');
ylabel('\Delta VC', 'FontSize', 17);
xlabel('Time/s')
savefig(gcf, [dir figName '.fig']);


