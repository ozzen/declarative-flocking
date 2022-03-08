function [num_violations, total_duration, res] = distanceViolations(trajs, params)
%distanceViolations - Computes the number of safety violations.
%
% Inputs: 
%    trajs - 1xn Array of structs containing the positions, velocites ...
%            for n runs.
%    params - struct containing the parameters.
%
% Outputs:
%    res - mx2 array of results. 
%      
% Author: Usama Mehmood, Graduate Student, Stony Brook University, NY, US
% email address: umehmood@cs.stonybrook.edu 
% Website: https://usamamehmood.weebly.com
% December 2018; Last revision: 05-Jan-2019
%------------- BEGIN CODE --------------
    res = [trajs.min_pw];
    res = find(res < params.dmin);
    violation_durations = diff(find(diff([nan ; res(:) ; nan]) ~= 1));
    
    num_violations = numel(violation_durations);
    total_duration = sum(violation_durations) * params.dt;
    res = [floor(res/params.steps)+1 mod(res, params.steps)];
end

