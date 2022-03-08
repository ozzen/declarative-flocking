function [ pos, vel ] = dynamics( pos, vel, acc, params )
%% dynamics() - dmpc
vel = vel + params.dt*acc;
for j = 1:size(pos,2)
    if params.predator && j == params.n
       if norm(vel(:,j)) > params.pFactor * params.vmax
            vel(:,j) = (params.pFactor * params.vmax / norm( vel(:,j) )) * vel(:,j);
       end
       continue;
    end
    if norm(vel(:,j)) > params.vmax
        vel(:,j) = (params.vmax/norm( vel(:,j) )) * vel(:,j);
    end
end
pos = pos + params.dt*vel;
end

