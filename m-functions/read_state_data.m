function [x, y, vx, vy] = read_state_data(file, numBirds)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Num = numBirds;
formatSpec = '%f';
sizeA = [4*Num Inf];
fileID = fopen(file, 'r');
A = fscanf(fileID,formatSpec,sizeA);
A = A';
fclose(fileID);
x = A(:, 1:Num);
y = A(:, Num+1:2*Num);
vx = A(:, 2*Num+1:3*Num);
vy = A(:, 3*Num+1:4*Num);
end

