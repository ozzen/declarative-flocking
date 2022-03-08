function [ res ] = CompositeBF(pos, vel, dmin, maxAcc )
    Num = size(pos, 2);
    res = 1;
    for i = 1:Num
        for j = i+1:Num
            pp = [pos(:,i)' ; pos(:,j)'];
            vv = [vel(:,i)' ; vel(:,j)'];
            res = res * barrierFunction( pp, vv, dmin, maxAcc );
        end
    end   
end

 