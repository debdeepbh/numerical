% get alpha from the original signal
function alpha = multi_getscOri_fpar(wax, aximp, testyori, type, p, noiseax)

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



% use bisection method
for j=1:p+1
	coeffvals = coeff(wori, p, j);
	PsiSq = abs(fft(B(j,:)));
	% initial guess
	a = 0.001;
	b = 1;
	while (abs(a-b)>errorl)
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

	alpha(j) = b;

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
			if i==k
				numera = numera + (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
			else
				denom = denom + (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
			end
		end
	end
	
	% sum over f
	Q(j) = sum(abs(fori).^4 .* PsiSq.^2 .* numera) / sum(abs(fori).^4 .* PsiSq.^2 .* denom);

	ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals)/ (1 + Q(j));
	%%%%%%%%%%% done calculating Q(j)
	% get the ratio of unthresholded coefficients

	fb = ratiounthr;

	
	% if the signs are opposite
	if ((a-fa)>=0 && (b - fb)<=0)
		% disect
		c = (a+b)/2;
		alpha(j) = c;

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
			if i==k
				numera = numera + (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
			else
				denom = denom + (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
			end
		end
	end
	
	% sum over f
	Q(j) = sum(abs(fori).^4 .* PsiSq.^2 .* numera) / sum(abs(fori).^4 .* PsiSq.^2 .* denom);

	%%%%%%%%%%% done calculating Q(j)
	% get the ratio of unthresholded coefficients
	ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals)/ (1 + Q(j));
		fc = ratiounthr;
		if ((c -fc) >0)
			a = c;
		else
			b = c;
		end
	% if the signs are the same
	elseif ((a-fa)<=0 && (b - fb)>=0)
		c = (a+b)/2;
		alpha(j) = c;

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
			if i==k
				numera = numera + (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
			else
				denom = denom + (noise(i)^2)*(noise(k)^2)*hsq(i,:)./(R(i,:).^2) ./ R(k,:);
			end
		end
	end
	
	% sum over f
	Q(j) = sum(abs(fori).^4 .* PsiSq.^2 .* numera) / sum(abs(fori).^4 .* PsiSq.^2 .* denom);

	%%%%%%%%%%% done calculating Q(j)
	ratiounthr = sum( coeffvals > sigmal)/ length(coeffvals)/ (1 + Q(j));
	% get the ratio of unthresholded coefficients
		fc = ratiounthr;
		if ((c -fc) >0)
			b = c;
		else
			a = c;
		end
	else
		c = 0.01;
		disp('No solution in between for');
		j
		break;
	end

	end
	alpha(j) = c;
end


