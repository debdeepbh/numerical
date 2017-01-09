function out = medF(nbdsize, f)
% median filter
% ex: medF(2,G)
f = double(f);
[m, n] = size(f);
out = zeros(size(m,n));

for i = 1:m
    for j = 1:n
       submat = f(max(1, i- nbdsize): min(m, i+ nbdsize), max(1, j- nbdsize): min(n, j+ nbdsize));
       out(i,j) = median(submat(:));
    end
end
out = uint8(out);