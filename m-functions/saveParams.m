function [ret] = saveParams(destination, params)
fid = fopen(destination, 'w');
fieldN = fieldnames(params);
fieldV = struct2cell(params);
for i = 1:numel(fieldN)
    fprintf(fid, "%-9s %-12f\n", [fieldN{i} ':'], fieldV{i});
end
ret = 1;
end

