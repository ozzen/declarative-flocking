N = GenMatricTwo(2, 2 );
n = size(N,1);
res = [];
for i = 1:n-1
    row = []
    for j = 1:n
        if j ~= 1
            row = [row '&' num2str(N(i, j)) ' '];
        else
            row = [row num2str(N(i, j)) ' '];
        end
    end
    row = [row ' \\'];
    res = [res ; row];
end

disp(res);



