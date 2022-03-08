function [ret] = murmuration_cost(pos, params)

n = size(pos,2);
ret = 0;
for i = 2:n
    pos_new = [pos(:,1) pos(:,i)];
    ret = ret + params.w_m(i-1) * (params.wc * average_squared_distance(pos_new))...
        + params.ws * separation(pos_new);
end
end

