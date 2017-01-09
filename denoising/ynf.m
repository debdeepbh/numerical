function out = ynf(nbdsize, variance, I)
% Works good with Gaussian noise instead of salt & pepper
% If the value of the noisy pixels are not too extreme
% example for salt&pepper: ynf(2,250,N); removes half of it
f = double(I);
[m, n] = size(f);
out = zeros(size(m,n));

for i = 1:m
    for j = 1:n
        s = 0;
        C = 0;
        for p = max(1, i- nbdsize):min(m, i+ nbdsize)
           for q =  max(1, j- nbdsize):min(n, j+ nbdsize)
               if abs(f(p,q) -f(i,j)) < variance
                   ker = 1;
               else
                   ker = 0;
               end
               s = s + ker*f(p,q);
               C = C + ker;
           end
        end
        out(i,j) = s/C;
          
    end
end
out = uint8(out);
