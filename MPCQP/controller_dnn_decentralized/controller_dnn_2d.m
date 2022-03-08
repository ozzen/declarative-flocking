function [a] = controller_dnn_2d(pos ,vel, i, params, net)
%% controller_dnn_2d - Run MPC controller and return first action
% Input:
%   - pos            % 2 x n - positions of the flock 
%   - vel            % 2 x n - velocities of the flock 
%   - i              % 1 x 1 - index of the central agent 
%   - params         % parameters
%                      params.knn should match the number of inputs of NN.\
%   - net            % Neural net
%
% Output:
%   - a              % 2 x 1 - The control action 
% Usama Mehmood - Dec 2019

a = zeros(2, 1);
% input_net = zeros(1, 4*(params.knn+1));

%% Load controller
% The caller loads it and passes it to this function.

%% Squared distances
sq_d = sq_distances_pairwise(pos);

% Find nearest neighbors
[~, Index] = sort(sq_d(i,:));
k_nearest_n = Index(1:params.knn+1);
% Re-arrange input to th NN - [Depends on architecture]
p = pos(:,k_nearest_n)';
p = p - p(1,:);
p = p(2:end,:);

v = vel(:,k_nearest_n)';
i_n = [p(:)', v(:)'];
% standard scaling
% input_net(i,:) = i_n;  

% scaling
% input_scaled = (input_net - mean(input_net)) ./ std(input_net,0,1);
% input_scaled = i_n;
a = predict(net,i_n)';
a = trim_vec(a, params.amax);

end

