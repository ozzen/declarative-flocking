function [] = plotDistanceBrake(posi, veli, velMagj, aMax, dt, dmin)
%plotDistanceBrake - Computes and plots the locus of points at braking
%distance apart from the ith Agent, for a jth agent that is moving towards
%it at a speed of velMagj.
%
% Inputs: 
%    posi - 1x2 Matric containing the 2D position.
%    veli - 1x2 Matric containing the 2D velocity.
%    velMagj - 1x3 The speed of the neighbor j. 
%    aMax - 1x1 Magnitude of maximum acceleration
%    dt - 1x1 The integration time-step for the simulation.
%    dmin - 1x1 Value for constraint on minimum pairwise distance.
%
% Outputs:
%      
% Author: Usama Mehmood, Graduate Student, Stony Brook University, NY, US
% email address: umehmood@cs.stonybrook.edu 
% Website: https://usamamehmood.weebly.com
% December 2018; Last revision: 05-Jan-2019
%------------- BEGIN CODE --------------
    theta = 0:1:360;
    temp = repmat(veli, numel(theta) ,1);
    temp2 = [cosd(theta)' sind(theta)'];
    vRel = sum( [cosd(theta)' sind(theta)'] .* repmat(veli, size(theta,1) ,1), 2 ) + velMagj;
    vRel(vRel<0) = 0;
    steps = floor((vRel) / (aMax * dt));
    
    d_br = dt * (steps .* vRel - 0.5 * (steps - 1) .* steps * aMax * dt);
%     d_br = vRel.^2 / (2 * aMax);
   
    d_br = d_br + dmin;
    resX = d_br.*cosd(theta') + posi(1);
    resY = d_br.*sind(theta') + posi(2);

    plot(resX, resY, "Color", [0, 1, 0], "LineWidth", 0.5);
%     plot(d_br_cont.*cosd(theta), d_br_cont.*sind(theta), "Color", [0, red, 0], "LineWidth", 1.5);
end

