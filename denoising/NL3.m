function out = NL3(nbdsize, variance, a, I) % Nonlocal Mean
% ex: NL3(5,20,3,G) without center mod

kize = 2;

f = double(I);
[m, n] = size(f);
out = zeros(size(m,n));

Ga = fspecial('gaussian',2*kize+1,a);
M = zeros(2*kize+1);

    function perExt = g(xval,yval)  % extension of f
        modxval = mod(xval,m);
        modyval = mod(yval,n);
        if modxval == 0
            modxval = m;
        end
        if modyval == 0;
            modyval = n;
        end
        perExt = f(modxval,modyval);
    end

for i = 1:m         %(i,j) = x
    txtlen = fprintf('Progress: %0.1f%%\n', i/m*100);
    for j = 1:n
        sSum = 0;
        C = 0;
        biggest = 0;
        for p = max(1, i- nbdsize):min(m, i+ nbdsize)   % (p,q) = y
            for q =  max(1, j- nbdsize):min(n, j+ nbdsize)
 %              if ~((p==i) ||(q==j))
                    % % can't express n2 = sum(sum(((g(-kize + p: kize + p: -kize + q : kize + q) - g(-kize + i: kize + i: -kize + j : kize + j)).^2) .* Ga));
                    for r = -kize:kize
                        for s = -kize:kize
                            M(1+kize+r,1+kize+s) = (g(p-r,q-s) - g(i-r,j-s))^2;
                        end
                    end
                    n2 = sum(sum(M .* Ga));
                   
                    
                    ker = exp( (-1)* n2 / 2/ (variance^2)) ;
%                    biggest = max(ker,biggest);
                    sSum = sSum + ker*f(p,q);
                    C = C + ker;

%                end
            end
        end

%         s = s + biggest*f(i,j);
%         C = C + biggest;
        out(i,j) = sSum/C;


    end
    fprintf(repmat('\b',1,txtlen));
end
fprintf('\n');
out = uint8(out);
end



