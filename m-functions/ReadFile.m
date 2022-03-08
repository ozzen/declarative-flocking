function [ ret] = ReadFile( name, format )
fileID = fopen(name,'r');
formatSpec = format;
ret = fscanf(fileID,formatSpec);
fclose(fileID);
end

