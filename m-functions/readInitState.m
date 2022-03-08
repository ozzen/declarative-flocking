function [x, y, vx, vy] = readInitState(init_path, simNum)
    p = sprintf(init_path, simNum);
    fileID = fopen(p,'r');
    formatSpec = '%f';
    sizeA = [1 Inf];
    A = fscanf(fileID,formatSpec,sizeA);
    n = A(1);
    A = reshape(A(2:end), [4, n])';
    x = A(:,1)';
    y = A(:,2)';
    vx = A(:,3)';
    vy = A(:,4)';
end

