function [bool] = FSL(pos, vel, params)
vel_rel = vel(:,1) - vel(:,2:end);
pos_rel = pos(:,1) - pos(:,2:end);

[h_pairwise, delta_v] = h_function(pos_rel, vel_rel, params);

%% Maximum change in bf possible in one time step.
num_steps = 1;
[switching_boundary] = max_change_boundary(pos_rel, vel_rel, params, num_steps);

%%
%     bool = any( h_pairwise(delta_v > 0) < switching_boundary);
bool = any( h_pairwise < switching_boundary);
end

