% get alpha from the original signal
% get the optimal scaling values for multi_fone()

function [alpha, ratioUnth, sigl] = multi_getscOri_fone(wax, aximp, testyori, type, p, sigma)

errorl = 0.0001;

N = length(wax);
hsq = abs(fft(aximp)).^2;
fori = fft(testyori);
fimp = fft(aximp);

wori = wtrans(testyori, type, p);
B = getbasismat(type, p, N);

alpha = zeros(1,p+1)+1; % start with scaling value 1
ratioUnth = zeros(1,p+1)+1; % start with scaling value 1
sigl = zeros(1,p+1)+1; % start with scaling value 1



% use bisection method
count = 1;
for i=1:p+1
	% get the wavelet coefficients of the level
	% of the original signal
	coeffvals = coeff(wori, p, i);
	PsiSq = abs(fft(B(i,:)));
	% initial guess
	a = 0.005;
	b = 0.5;
	while (abs(a-b)>errorl)
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

		alpha(i) = b;
		%%%%%%%%%%%
		% using alpha and sigma, get multiplier lambda
		mult = hsq ./( hsq + N*alpha(i)*(sigma^2)./(abs(fori).^2) );
		% get level dependent sigmal
		sigmal = sqrt( dot( (PsiSq./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);
		% get the ratio of original wavelet coefficients that are bigger than the noise value
		ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals);
		%%%%%%%%%%5
		fb = ratiounthr;

		
		% if the signs are opposite
		if ((a-fa)>=0 && (b - fb)<=0)
			% disect
			c = (a+b)/2;
			alpha(i) = c;
		%%%%%%%%%%%
		% using alpha and sigma, get multiplier lambda
		mult = hsq ./( hsq + N*alpha(i)*(sigma^2)./(abs(fori).^2) );
		% get level dependent sigmal
		sigmal = sqrt( dot( (PsiSq./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);
		% get the ratio of unthresholded coefficients
		ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals);
		%%%%%%%%%%5
			fc = ratiounthr;
			if ((c -fc) >0)
				a = c;
			else
				b = c;
			end
		% if the signs are the same
		elseif ((a-fa)<=0 && (b - fb)>=0)
			c = (a+b)/2;
			alpha(i) = c;
		%%%%%%%%%%%
		% using alpha and sigma, get multiplier lambda
		mult = hsq ./( hsq + N*alpha(i)*(sigma^2)./(abs(fori).^2) );
		% get level dependent sigmal
		sigmal = sqrt( dot( (PsiSq./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);
		% get the ratio of unthresholded coefficients
		ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals);
		%%%%%%%%%%5
			fc = ratiounthr;
			if ((c -fc) >0)
				b = c;
			else
				a = c;
			end
		else
			%c = 0.001;
			c = 1;
			disp('No solution in between for');
			i
			break;
		end
	count = count + 1;
	end
	alpha(i) = c;
	ratioUnth(i) = ratiounthr;
	sigl(i) = sigmal;
	
end


