clc
clear all
close all

%% Create and save a rectangular obstacles in a text file
% Output file format:
%
%   number of obstacles
%   centre x-coordinate for first rectangle
%   centre y-coordinate for first rectangle
%   width for first rectangle
%   length for first rectangle
%   angle of orientation, set to zero
%   .
%   .
%   .
%
% Method:
%   1. Draw n rectangles on the figure window.
%   2. Then press enter key.
%
%% Inputs
n = 10;
target = [60, 50];
output_file = 'rectangles1.txt';

%%
figure;

plot(target(1), target(2), 'k.', 'MarkerSize', 20);
axis equal
axis([-15 100 -15 100]);
rectangle('position', [-15 -15 30 30]);


for i = 1:n
    h(i) = imrect(); %Takes input a rectangle drawn on the figure window
end

currkey=0;
% do not move on until enter key is pressed
while currkey~=1
    pause; % wait for a keypress
    currkey=get(gcf,'CurrentKey'); 
    if currkey=='return'
        currkey=1;
    else
        currkey=0;
    end
end

pos = zeros(n, 4);
for i = 1:n
    pos(i,:) = h(i).getPosition;
end

centre = pos(:,1:2) + 0.5*pos(:,3:4);
wl = pos(:,3:4);

fileID = fopen(output_file, 'w');
fprintf(fileID, '%d\n', n);
for i =  1:n
    temp =  [centre(i,:) wl(i,:) 0];
    fprintf(fileID, '%f\n%f\n%f\n%f\n%f\n', temp');
end