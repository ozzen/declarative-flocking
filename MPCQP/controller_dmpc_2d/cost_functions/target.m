function [res] = target(pos, target)
%% dmpc
res =  norm(pos(:,1) - target);
end

