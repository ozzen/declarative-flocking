function [s,x,y,z,vx,vy,vz,accx,accy,accz,omegas,params] = read_data_mat(traj, is_quad)
s = traj.s;
x = traj.x;
y = traj.y;
z = traj.z;
vx = traj.vx;
vy = traj.vy;
vz = traj.vz;
accx = traj.accx;
accy = traj.accy;
accz = traj.accz;
params = traj.params;
if params.quad
    omegas = traj.omegas;
else
    omegas = 0;
end
end

