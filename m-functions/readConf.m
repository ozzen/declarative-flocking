function [params ] = readConf( fname )
fid = fopen(fname,'r');
while ~feof(fid)
    tline = fgets(fid);
    Data = strsplit(tline);
    if Data{1} == "amax"
        params.amax = str2double(Data{2});
        continue;
    end
    if Data{1} == "vmax"
        params.vmax = str2double(Data{2});
        continue;
    end
    if Data{1} == "dt"
        params.dt = str2double(Data{2});
        continue;
    end
    if Data{1} == "ct"
        params.ct = str2double(Data{2});
        continue;
    end
    if Data{1} == "dmin"
        params.dmin = str2double(Data{2});
        continue;
    end
    if Data{1} == "alpha"
        params.alpha = str2double(Data{2});
        continue;
    end
    if Data{1} == "learningRate"
        params.lr = str2double(Data{2});
        continue;
    end
    if Data{1} == "maxIters"
        params.maxiters = str2double(Data{2});
        continue;
    end
    if Data{1} == "horizon"
        params.h = str2double(Data{2});
        continue;
    end
    if Data{1} == "steps"
        params.steps = str2double(Data{2});
        continue;
    end
    if Data{1} == "rc"
        params.rad_sensing = str2double(Data{2});
        continue;
    end

end
end

