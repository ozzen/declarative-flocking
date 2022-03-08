% from Olfati-saber, but it's different from J_v of Zhang2015
function K = velocityMismatch(p)
K=sum(sum(p.^2,2))/2/size(p,1);
end