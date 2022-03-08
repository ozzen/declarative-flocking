clc
clear all
%%
figure;
n = 25;

obstBB = [60 -30 100 60]; %[x, y, l, w]
areaObstBB = obstBB(3)*obstBB(4);

points = [obstBB(3)*rand(n, 1) + obstBB(1), obstBB(4)*rand(n, 1) + obstBB(2)];

% points(end,:) = [180, 0];

rectangle('position', [-15 -15 30 30]);
axis equal
axis([-30 160 -30 30])
hold on
plot(points(:,1),points(:,2), '.' ,'MarkerSize', 5)

% points = zeros(n,2);
% for i = 1:n
%     if i == n
%         title('target')
%         [x, y] = ginput(1);
%         plot(x,y,'r*', 'MarkerSize', 6);
%         points(i,1) = x;
%         points(i,2) = y;
%         continue;
%     end
%     [x, y] = ginput(1);
%     points(i,1) = x;
%     points(i,2) = y;
%     plot(x,y,'b*', 'MarkerSize', 3);
% end

%% Save
fileID = fopen(['pointObstacles' num2str(n) '.txt'],'w');
fprintf(fileID, '%d\n', n);
for i =  1:n
    fprintf(fileID, '%f\n%f\n', points(i,:));
end
