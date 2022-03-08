function [a] = controller_dnn_3d(s, params, net)
%% controller_dnn_3d - Run MPC controller and return first action
% Non-liear constraints passed to fmincon as function handle. 
% like this: 
%     x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon)
% Input:
%   - s            % 1x12xn - state of the flock 
%   - params      d % parameters
%
% Output:
%   - a            % 3 x n - The sequence of control action over the horizon 
% Usama Mehmood - Oct 2019

a = zeros(3, params.n);
input_net = zeros(params.n, 6*(params.knn+1));

%% find neighbors and arrange inputs to the controller.
pos = squeeze(s(1,1:3,:));
vel = squeeze(s(1,4:6,:));

sq_d = sq_distances_pairwise(pos);

for i = 1:params.n
% Find nearest neighbors
    [~, Index] = sort(sq_d(i,:));
    k_nearest_n = Index(1:params.knn+1);
% Re-arrange input to th NN - [Depends on architecture]
    p = pos(:,k_nearest_n)';
    v = vel(:,k_nearest_n)';
    i_n = [p(:)', v(:)'];
% standard scaling
    input_net(i,:) = i_n;  
end

% scaling
% input_scaled = (input_net - mean(input_net)) ./ std(input_net,0,1);
input_scaled = input_net;
%% Run controller for all agents
for i = 1:params.n
    a(:,i) = predict(net,input_scaled(i,:))';
    a(:,i) = trim_vec(a(:,i), params.amax);
end

end

