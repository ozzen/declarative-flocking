function [min_distance, closest_point, count] = minFlockObstacleDistance(q, obs)
    min_distance = inf;
    n = size(q, 1);
    closest_point = zeros(n, 2);
    count = 0;
    for k = 1:n
        [dist, closest_point(k,:)] = point2rects(q(k,:), obs);
        if dist < min_distance
            min_distance = dist;
        end
        if dist < 1
            count = count + 1;
        end
    end
end

