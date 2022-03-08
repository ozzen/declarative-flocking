function [Rect] = readRects(filename)

fileID = fopen(filename,'r');
ret = fscanf(fileID,"%f");
n = ret(1);
rdata = reshape( ret(2:end) , [5, n] )';

%% Add to struct.
for i = 1:n
    Rect(i).cx =rdata(i, 1);
    Rect(i).cy = rdata(i, 2);
    Rect(i).w = rdata(i, 3);
    Rect(i).l = rdata(i, 4);
    Rect(i).theta = rdata(i, 5);
end

end

