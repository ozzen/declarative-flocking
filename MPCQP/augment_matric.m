function [ret] = augment_matric(M)
s = size(M);

A = zeros(s(1)^2, s(2));
for i = 1:s(1)
    A(s(1) * i - s(1) + 1:s(1) * i,:) = repmat(M(i,:), s(1), 1);
end
 

B = repmat(M, s(1), 1);
ret = [A, B];

end

