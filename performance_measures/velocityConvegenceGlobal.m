function [ret] = velocityConvegenceGlobal(p, params)
v_bar = mean(p);
% v_dev = 0;
% for i = 1:params.N
%     v_dev = v_dev + norm(p(i,:) - v_bar);
% end

ret = sum(sum( (p - v_bar).^2, 2))/params.N;

% ret = v_dev / params.N;
end

