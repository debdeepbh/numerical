function out = gaussF2(nbdsize, variance, I)    % Gaussian filter
% ex: gaussF2(5,2,G)
f = double(I);
[m, n] = size(f);
out = zeros(size(m,n));

for i = 1:m
    for j = 1:n
        s = 0;
        C = 0;
        for p = max(1, i- nbdsize):min(m, i+ nbdsize)
           for q =  max(1, j- nbdsize):min(n, j+ nbdsize)
               n2 = (p-i)^2 + (q-j)^2;
               ker = exp( (-1)* n2 / 2/ (variance^2)) ;
               s = s + ker*f(p,q);
               C = C + ker;
           end
        end
        out(i,j) = s/C;
          
    end
end
out = uint8(out);
