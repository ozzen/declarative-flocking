function [amax, vmax, beta] = ParameterRanges(source)
    listing = dir(source);
    folders = extractfield(listing, 'name');
    
    num = numel(folders) - 3;
    amax = zeros(1, num);
    vmax = zeros(1, num);
    beta = zeros(1, num);

    for i = 1:num-1
        experiment_name = folders{i + 3}
        toks = split(experiment_name, '_');
        amax(i) = str2double(toks{2});
        vmax(i) = str2double(toks{3});
        beta(i) = str2double(toks{4});
    end
    amax = unique(amax);
    vmax = unique(vmax);
    beta = unique(beta);
    
    amax = amax(2:end);
    vmax = vmax(2:end);
    beta = beta(2:end);

end

