function alpha = multi_plotalpha(wax, aximp, testyori, type, p, sigma)

% for j = 1:p+1
for j = 1:p+1;

N = length(wax);
hsq = abs(fft(aximp)).^2;
fori = fft(testyori);
fimp = fft(aximp);

B = getbasismat(type, p, N);
PsiSq = abs(fft(B(j,:)));

i=1;
for alpha=0.01:0.05:1
	% using alpha and sigma, get multiplier lambda
	mult = hsq ./( hsq + N*alpha*(sigma^2)./(abs(fori).^2) );
	% get level dependent sigmal
	sigmal = sqrt( dot( (PsiSq./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);

	% get the ratio of unthresholded coefficients
	coeffvals = coeff(wtrans(testyori, type, p), p, j);
	ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals);

	xval(i) = alpha;
	yval(i) = alpha - ratiounthr;
	i = i+1;
end


plot(xval,yval)
hold on
end

% the right hand side
