function [res] = vm(s, params)
vel = squeeze(s(1,4:6,:));
rel_speed = sqrt(sum((repmat(vel(:,1),1,params.knn) - vel(:,2:end)).^2, 1));

%% AWNing
res = sum(params.w_m(1:params.knn) .* rel_speed);
end

