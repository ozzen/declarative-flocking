function [ret] = augment_matric_pair(M, N)
s_m = size(M);
s_n = size(N);

B = repmat(N, s_m(1), 1);

A = zeros(s_m(1) * s_n(1), s_m(2));
for i = 1:s_m(1)
    A(s_n(1) * i - s_n(1) + 1:s_n(1) * i,:) = repmat(M(i,:), s_n(1), 1);
end

ret = [A, B];

end

