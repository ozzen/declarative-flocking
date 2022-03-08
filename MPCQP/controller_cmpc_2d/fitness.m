function [res] = fitness(pos, vel, params)
%#codegen
%% cmpc - fitness 
res = 0;
%% Obstacle Avoidance
% rects = params.rects;
% obst = obstacle_cost(pos, rects, params);

%% Target seeking
% tar = [60; 50];
% tgt = target(pos, tar, params);

%% separation
sp = separation_polynomial(pos);
% spC = separation_constrained(pos, params);

%% cohesion
asd = sum_sq_distances( pos );

%% vm
% vm = velocity_matching(vel, params) ;

%% Predator avoiadance
if params.predator
    [res] = predator_avoidance(pos, params);
end

% res = asd + 1000000 * spC;
res = params.wc * asd...
    + params.ws * sp;
%     + params.wo * obst...
%     + params.wt * tgt;

end

% tgt weight: 1, 10, 50 