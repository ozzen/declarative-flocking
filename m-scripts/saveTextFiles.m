function [] = saveTextFiles(sourceRes, sourceInit, destination, confFile, simNum)

initFile = sprintf(sourceInit, simNum);

x_file = sprintf(sourceRes, simNum, 'x');
y_file = sprintf(sourceRes, simNum, 'y');
vx_file = sprintf(sourceRes, simNum, 'vx');
vy_file = sprintf(sourceRes, simNum, 'vy');
ax_file = sprintf(sourceRes, simNum, 'ax');
ay_file = sprintf(sourceRes, simNum, 'ay');
f_file = sprintf(sourceRes, simNum, 'fitness');

mkdir(destination);

copyfile(initFile, destination);

copyfile(x_file, destination);
copyfile(y_file, destination);
copyfile(vx_file, destination);
copyfile(vy_file, destination);
copyfile(ax_file, destination);
copyfile(ay_file, destination);
% copyfile(f_file, destination);
copyfile(confFile, destination);

end