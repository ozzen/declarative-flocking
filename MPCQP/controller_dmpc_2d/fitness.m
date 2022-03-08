function [res] = fitness(pos, params)
%#codegen
%% dmpc - fitness
% asd = average_squared_distance(pos);
% sp = separation(pos);
%% Obstacle avoidance
% rects = params.rects;
% obst = obstacle_cost(pos, rects, params);

%% Target seeking
% tar = [60; 50];
% tgt = target(pos, tar);

%% Predator Avoidance
% res = predator_avoidance(pos, params);

% weighted sum
% res = params.wc * asd ...
%     + params.ws * sp;
%     + params.wo * obst ...
%     + params.wt * tgt;

%% murmuration
res = murmuration_cost(pos, params);
end

