function [ pos, vel ] = Dynamics( pos, vel, acc, params )
    vel = vel + params.dt*acc;
%     Map = sqrt( vel(1,:).^2 + vel(2,:).^2 ) > params.vmax * ones(1, params.n);
%     temp = sqrt( vel(1,:).^2 + vel(2,:).^2 ); 
%     temp = params.vmax * ones(size(temp)) ./ temp;
%     vel = vel .* [Map; Map] .* [temp; temp] + vel .* [~Map; ~Map];
    for j = 1:params.n
        if norm(vel(:,j)) > params.vmax
            vel(:,j) = (params.vmax/norm( vel(:,j) )) * vel(:,j);
        end
    end
    pos = pos + params.dt*vel;
end

