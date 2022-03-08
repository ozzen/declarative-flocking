function [str] = num2string(num)
%num2string - Outputs string represintation of number correct to 2 decimal places.
%
% Inputs: 
%    num - 1x1 double.
%
% Outputs:
%    str - 1x1 string.
%    examples: 1->1.00
%              1.1->1.10 
%              1.15->1.15 

%      
% Author: Usama Mehmood, Graduate Student, Stony Brook University, NY, US
% email address: umehmood@cs.stonybrook.edu 
% Website: https://usamamehmood.weebly.com
% December 2018; Last revision: 05-Jan-2019
%------------- BEGIN CODE --------------
if any(num == 1:9)
    str = [num2str(num) '.00'];
elseif mod(uint16(num * 100),  10) < 10e-6
    str = [num2str(num) '0'];
else
    str = num2str(num);
end




end

