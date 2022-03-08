clc
clear all;
%%
vMag = 1.0;
aMax = 1.625;
dt = 0.10;
dmin = 3;

% isoVs = 0:0.5:2;
% theta = 0:1:360;
% red = 0;
% for isoV = isoVs
%     vRel = vMag * cosd(theta) + isoV;
%     vRel(vRel<0) = 0;
%     steps = floor((vRel) / (aMax * dt));
%     d_br = dt * (steps .* vRel - 0.5 * (steps - 1) .* steps * aMax * dt);
%     
%     d_br_cont = vRel.^2 / (2 * aMax);
%     d_br_cont = d_br_cont + dmin;
%     
%     d_br = d_br + dmin;
%     red = red + 1/(numel(isoVs)+1);
%     plot(d_br.*cosd(theta), d_br.*sind(theta), "Color", [red, 0, 0], "LineWidth", 1.5);
%     plot(d_br_cont.*cosd(theta), d_br_cont.*sind(theta), "Color", [0, red, 0], "LineWidth", 1.5);
% 
%     hold on
%     axis equal
% end

plotDistanceBrake([2,2], [-1,1], vMag, aMax , dt, dmin);
axis equal;

