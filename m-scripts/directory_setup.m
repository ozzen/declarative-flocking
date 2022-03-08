% Setup the following:
% 1. Source path
% 2. Destination path
% 3. Results folder

%% NIPS 2019
simNum = 2;
src_folder = 'exp7';
dest_folder = 'src1';
g_d = 'g';

root_dir = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/';
prefix = ['experiments_small_flock/' src_folder '/'] ; %Source path
% prefix = '/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/';

destPath = [root_dir prefix 'Results/'];
sourcePath = [root_dir prefix g_d 'mpc_%d_%s.txt'];
mkdir(destPath);

initPathName = [root_dir 'experiments_small_flock/' dest_folder '/init_conf_%d.txt'];
confFile = [root_dir prefix 'result_files/' g_d 'mpc.conf'];
% obst_file = [root_dir '/rectangles.txt'];
obst_file = [root_dir '/results_neurips/obstacles/1/rectangles.txt'];

traj_file = [root_dir prefix 'result_files/traj_data.mat'];
