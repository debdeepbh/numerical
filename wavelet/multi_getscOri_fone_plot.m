% plot ratiounthr vs alpha_j
% get alpha from the original signal
% get the optimal scaling values for multi_fone()

function alphaout = multi_getscOri_fone_plot(wax, aximp, testyori, type, p, sigma)

errorl = 0.01;

N = length(wax);
hsq = abs(fft(aximp)).^2;
fori = fft(testyori);
fimp = fft(aximp);

wori = wtrans(testyori, type, p);
B = getbasismat(type, p, N);

alpha = zeros(1,p+1)+1; % start with scaling value 1
ratioUnth = zeros(1,p+1)+1; % start with scaling value 1
sigl = zeros(1,p+1)+1; % start with scaling value 1


tic
% use bisection method
for i=1:p+1
	% get the wavelet coefficients of the level
	% of the original signal
	coeffvals = coeff(wori, p, i);
	PsiSq = abs(fft(B(i,:)));
	% initial guess
	%a = 0.0001;

	%minval(i) = 1;
	count = 1;
	for a=1e-3:1e-3:1	
		alphaval(count) = a;
		
		alpha(i) = a;
		%%%%%%%%%%%
		% using alpha and sigma, get multiplier lambda
		mult = hsq ./( hsq + N*alpha(i)*(sigma^2)./(abs(fori).^2) );
		% get level dependent sigmal
		sigmal = sqrt( dot( (PsiSq./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);
		% get the ratio of unthresholded coefficients
		ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals);
		%%%%%%%%%%5
		fa = ratiounthr;

		ratval(count) = fa;

		diffvec(count) = abs(a-fa);
		notdiffvec(count) = a - fa;


		count = count + 1;
	end
	%figure 1;
	%subplot(p+1,1,i);
	%plot(alphaval, log(alphaval), alphaval, log(ratval));
	%legend('alpha', 'f(alpha)');

	%figure 2;
	%subplot(p+1,1,i);
	%plot(alphaval, log(diffvec));


	%sigl(i) = sigmal;
	
	alphaout(i) = alphaval(find(diffvec == min(diffvec)))(1);

	
end


toc

