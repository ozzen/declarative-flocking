function [x, y, vx, vy, ax, ay, f, policy] = read_data(sim_number, num_birds, file_template)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x_file = sprintf(file_template, sim_number, 'x');
y_file = sprintf(file_template, sim_number, 'y');
vx_file = sprintf(file_template, sim_number, 'vx');
vy_file = sprintf(file_template, sim_number, 'vy');
ax_file = sprintf(file_template, sim_number, 'ax');
ay_file = sprintf(file_template, sim_number, 'ay');
f_file = sprintf(file_template, sim_number, 'fitness');
policy_file = sprintf(file_template, sim_number, 'policy');


x = read_data_internal(x_file, num_birds);
y = read_data_internal(y_file, num_birds);
vx = read_data_internal(vx_file, num_birds);
vy = read_data_internal(vy_file, num_birds);
ax = read_data_internal(ax_file, num_birds);
ay = read_data_internal(ay_file, num_birds);
f = read_data_internal(f_file, 1);
policy = read_data_internal(policy_file, 1);

    function [z] = read_data_internal(filename, cols)
        fp = fopen(filename, 'r');
        if fp > 0
            z = fscanf(fp, '%f', [cols Inf]);
            z = z';
            fclose(fp);
        else
            z = [];
        end
    end

end

