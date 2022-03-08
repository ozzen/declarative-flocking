function [] = plotPolicyAcc(accx, accy, policy, i)
%plotPolicyAcc - Plot of magnitude of acceleration and the control policy.
%
% Inputs: 
%    accx - 1xsteps double.
%    accy - 1xsteps double.
%    policy - 1xsteps double. (0/1)
%    i - 1x1 double.

%
% Outputs:
%      
% Author: Usama Mehmood, Graduate Student, Stony Brook University, NY, US
% email address: umehmood@cs.stonybrook.edu 
% Website: https://usamamehmood.weebly.com
% December 2018; Last revision: 05-Jan-2019
%------------- BEGIN CODE --------------
figure;
xLabel = 'Time';
yLabel = 'Acceleration and Policy';
Num = size(accx, 2);

A = sqrt(accx(2:end,i).^2 + accy(2:end,i).^2);
plotCurve(A, xLabel, yLabel, 1.5, Num);
hold on
plotCurve(policy(:,i), xLabel, yLabel, 1.5, Num)
end

