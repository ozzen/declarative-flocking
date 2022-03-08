clc
clear all
close all

addpath('utils')
addpath('performance_measures');
addpath('m-functions');
addpath('m-scripts');
addpath('Plotting')

%% Params
[params ] = readConf( "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/C Code/dmpc.conf" );
params.Tfin = (params.steps - 1) * params.dt;
params.N = 30; % number of agents
params.nd = 2; % space dimension, 2 = 2-D, 3 = 3-D

%%
source = 'experiments_loop/';
[amax, vmax, beta] = ParameterRanges(source);
numberOfViolations = zeros(numel(amax), numel(vmax), numel(beta));
durationOfViolations = zeros(numel(amax), numel(vmax), numel(beta));

total_runs = numel(amax) * numel(vmax) * numel(beta); 
count = 1;
for i = 1: numel(amax)
    for j = 1: numel(vmax)
        for k = 1: numel(beta)
            
            experiment_name = ['exp_' num2string(amax(i)) '_' num2string(vmax(j)) '_' num2string(beta(k))];
            %% Read Trajs
            nReps = 1 ;
            dataFolder = [source experiment_name '/'];
            str1 = 'dmpc_%d_%s.txt';
            stateFile = [dataFolder, str1];
            disp([num2str(count) ' out of ' num2str(total_runs)]);
            [ trajs] = makeTrajs( params , nReps, stateFile);
            
            [numOfViol, durOfViol] = distanceViolations(trajs, params);
            numberOfViolations(i,j,k) = numOfViol;
            durationOfViolations(i,j,k) = durOfViol;
            count = count + 1;
        end
    end
end

%% Plotting
total_plots = numel(beta);
[x, y] = meshgrid(amax, vmax);
map = [fliplr(linspace(0,0.8,100))' linspace(0,0.8,100)' zeros(100,1)];


for i = 0:total_plots - 1 %numel(beta)-1
    subplot(total_plots, 2, 2*i+1);
    surf(x, y, flip(flip(numberOfViolations(:,:,i+1), 2)', 1) );
    xlabel('amax');
    ylabel('vmax');
    zlabel('violations');
    title(['beta = ' num2string(beta(i+1))])
    view(0, 90);
    cb = colorbar;
    cb.Label.String = 'number of violations';

    
    subplot(total_plots, 2, 2*i+2);
    surf(x, y, flip(flip(durationOfViolations(:,:,i+1),2)', 1) );
    xlabel('amax');
    ylabel('vmax');
    zlabel('duration');
    title(['beta = ' num2string(beta(i+1))])
    view(0,90)
    cb = colorbar;
    cb.Label.String = 'duration / seconds';
end

colormap summer;
% rs = 2 * (2 * y).^2  ./ (4 * x);
% surf(x, y, rs)
figure;
plot_parameters(numberOfViolations, vmax, amax, beta)