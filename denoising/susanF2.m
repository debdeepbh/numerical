function out = susanF2(nbdsize, variance, I)
% works good with gaussian noise instead of salt&pepper
% example: susanF2(5,60,G)
f = double(I);
[m, n] = size(f);
out = zeros(size(m,n));

for i = 1:m
    for j = 1:n
       submat = f(max(1, i- nbdsize): min(m, i+ nbdsize), max(1, j- nbdsize): min(n, j+ nbdsize));
       ker = exp( -((submat - f(i,j)).^2) ./((variance)^2));
       out(i,j) = sum(sum( ker .* submat)) / sum(sum(ker));
    end
end
out = uint8(out);
