function [init] = gen_init(params, numSim)
% Usama Mehmood - Oct 2019
    
    init = cell(1,numSim);
    i = 1;
    while i <= numSim
        pos = params.ipos(1) + (params.ipos(2)-params.ipos(1)) * rand(1,3, params.n);
        vel = params.ivel(1) + (params.ivel(2)-params.ivel(1)) * rand(1,3, params.n);
        rest = zeros(1,6,params.n);
        s_init = cat(2, pos, vel, rest);
        if safe_dist_check(pos, params)
            init{i} = s_init;
            i = i + 1;
        end           
    end
end

