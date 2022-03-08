function [posN, velN, neighbors] = get_neighbors_states(ai, pos, vel, radius)
    distance_to_ai = sqrt(sum((pos - pos(:,ai)).^2, 1));
    neighbors = find(distance_to_ai < radius & distance_to_ai > 0 );

    posN = [pos(:,ai), pos(:, neighbors)];
    velN = [vel(:,ai), vel(:, neighbors)];
end

