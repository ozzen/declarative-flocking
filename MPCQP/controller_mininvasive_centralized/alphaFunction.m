function [ res ] = alphaFunction( CBF, gamma )
%     if CBF >= 100
%         res = gamma * (CBF.^3);
%     else
%         res = -1.0;
%     end  
    res = gamma * (CBF.^2);
% %     res = 0;
% res = -3;
end

