% get alpha from the original signal
function alphaout = multi_getscOri_fpar_plot(wax, aximp, testyori, type, p, noiseax)

errorl = 0.01;

[M, N] = size(wax);
% fft of a martrix is calculated column wise
hsq = (abs(fft(aximp')).^2)';
fori = fft(testyori);
fimp = fft(aximp')';
noise = std(noiseax')';

wori = wtrans(testyori, type, p);
B = getbasismat(type, p, N);

alpha = zeros(1,p+1)+1; % start with scaling value 1
ratioUnth = zeros(1,p+1)+1; % start with scaling value 1
sigl = zeros(1,p+1)+1; % start with scaling value 1



% use bisection method
for j=1:p+1
	coeffvals = coeff(wori, p, j);
	PsiSq = abs(fft(B(j,:)));
	% initial guess
	count = 1;
	for a = 1e-3:1e-3:1
		alphaval(count) = a;

		alpha(j) = a;
		%%%%%%%%%%%%%55 this is level dependent sigmal%%%%%%%%%%%55%%55
		for i=1:M
			mult = hsq(i,:) ./( hsq(i,:) + N*alpha(j)*(noise(i)^2)./(abs(fori).^2) );
			sigmalall(i) = sqrt( dot( (PsiSq./abs(fimp(i,:))).^2, abs(mult).^2 * noise(i)^2 )/N);
		end
		sigmal = sqrt(sum((sigmalall).^2))/M;	% avg over i
		%%%%%%%%%%%%%%%  calculate Q(j)
		% define Rij
		for i =1:M
			R(i,:) = abs(fori).*hsq(i,:) + N*noise(i)^2*alpha(j);
		end
		% inner sum
		denom = zeros(1,N);
		numera = zeros(1,N);
		for k=1:M
			for i=1:M
				sumterm = (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
				if i==k
					numera = numera + sumterm;
				else
					denom = denom + sumterm;
				end
			end
		end
		
		% sum over f
		Q(j) = sum(abs(fori).^4 .* PsiSq.^2 .* numera) / sum(abs(fori).^4 .* PsiSq.^2 .* denom);

		%%%%%%%%%%% done calculating Q(j)
		% get the ratio of unthresholded coefficients
		ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals)/ (1 + Q(j));
		%%%%%%%%%%5
		fa = ratiounthr;


		ratval(count) = fa;

		diffvec(count) = abs(a-fa);
		notdiffvec(count) = a - fa;


		count = count + 1;

	end
	%figure 1;
	%subplot(p+1,1,j);
	%plot(alphaval, log(alphaval), alphaval, log(ratval));
	%legend('alpha', 'f(alpha)');

	%figure 2;
	%subplot(p+1,1,j);
	%plot(alphaval, log(diffvec));


	alphaout(j) = alphaval(find(diffvec == min(diffvec)))(1);
end


