function [res] = predator_avoidance(pos, params)
%#codegen
%% dmpc - predator_avoidance

%% cohesion and separation of the flock
sp = separation(pos(:,1:end-1));
asd = average_squared_distance(pos(:,1:end-1));

%% predator cost

sq_dist_to_pred = sqrt(sum((pos(:,1) - pos(:,end)).^2, 1)).^3;
pa =  1/sq_dist_to_pred;

res = params.wc * asd...
    + params.ws * sp...
    + params.wp * pa;

end

