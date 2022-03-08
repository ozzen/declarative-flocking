function [ res ] = barrierFunction( pos, vel, dmin, maxAcc )
p_ij = pos(1,:) - pos(2,:);
v_ij = vel(1,:) - vel(2,:);
ret = sum((p_ij/norm(p_ij)).*(v_ij)) + sqrt( 4 * maxAcc * (norm(p_ij) - dmin) );
res = max(ret, 0);
% if norm(p_ij) > 4*dmin
%     res = 1;
% end

end

