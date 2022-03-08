% Gather statistics for centralized declarative flocking controller

addpath('../utils')
addpath('../Zhang_distributed_MPC/performance_measures')
clear;
params.nd = 2; % space dimension, 2 = 2-D, 3 = 3-D
params.N = 30; % number of agents
params.r = 8.4;
params.dt = 0.3;
params.Tfin = 100;
params.steps = 334;
params.d = 7;

trajs = [];
nReps = 100;
%dataFileTemplate = 'experiments/exp1/reynolds_%d_%s.txt';
dataFileTemplate = 'experiments/exp1/gmpc_%d_%s.txt';
%dataFileTemplate = 'experiments/exp3/gmpc_%d_%s.txt';
%dataFileTemplate = 'experiments/exp4-dmpcd/dmpcd_%d_%s.txt';
%dataFileTemplate = 'experiments/exp3/gmpcd_%d_%s.txt';
for rept = 1:nReps
    [x, y, vx, vy] = read_data(rept, params.N, dataFileTemplate);
    if params.steps < size(x, 1)
        x = x(1:params.steps,:);
        y = y(1:params.steps,:);
        vx = vx(1:params.steps,:);
        vy = vy(1:params.steps,:);
    end
    traj.q = zeros(params.steps,params.N,params.nd);
    traj.p = zeros(params.steps,params.N,params.nd);
    traj.pnet = zeros(params.steps,params.N,params.N);
    traj.connectivity = zeros(params.steps, 1);
    traj.v_converg = zeros(params.steps, 1);
    traj.irreg = zeros(params.steps, 1);
    traj.diameter = zeros(params.steps, 1);
    for t = 1:params.steps
        q = [x(t,:)' y(t,:)'];
        p = [vx(t,:)' vy(t,:)'];
        traj.q(t,:,:) = q;
        traj.p(t,:,:) = p;
        pnet = proximityNet(q, params.r);
        traj.pnet(t,:,:) = pnet;
        traj.connectivity(t) = connectedComponents(pnet);
        traj.v_converg(t) = velocityConvergence(p);
        %traj.irreg(t) = irregularity(q, params.r, params.d, pnet);
        traj.irreg(t) = newIrregularity(q, pnet);
        traj.diameter(t) = componentWiseDiameter(q, pnet);
    end
    %traj.diameter = flock_diameter(x, y, []);
    %traj.dist = sum(cumulative_travelled_distance(x, y, []),2);
    
    trajs = [trajs,traj];
end
plot_measures;