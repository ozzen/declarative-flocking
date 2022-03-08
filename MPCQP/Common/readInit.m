function [ posi, veli ] = readInit(file_template, sim_number )
%
    fname = sprintf(file_template, sim_number);
    fp = fopen(fname, 'r');
    
    N = str2double(fgetl(fp));
    A = fscanf(fp,'%f',[4 N]);
    posi = A(1:2,:);
    veli = A(3:4,:);
end

