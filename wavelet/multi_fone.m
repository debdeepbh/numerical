function [z, ratiounthres] = multi_fone(wax, aximp, testyori, type, p, sigma, sc_ax, rho, method)



[M, N] = size(wax);
% get the wavelet basis
B = getbasismat(type, p, N);


% all the antennas
%% get the  avg of the antennas, wtrans 
	fsig = fft(wax);
	fimp = fft(aximp);
	fori = fft(testyori);
	% to be chosen later
	scaling = sc_ax;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% estimate the wavelet coefficients 
	%% and the noise variance of summed data

	% matlab allows empty matrix []
	w = [];

	% varl stores the level dependent leaked noise variance

	% for level j
	for j=1:p
		% do deconvolution using specific scaling for the level
		% for vector valued noise
		
		[fdec, mult] = multi_fdecwien(fsig, fimp, sigma, fori, scaling(j));

		%[fdec, mult] = fdecwien(fsig, fimp, sigma_scalar, scaling(j));

		% prepare beta(k), the estimate for the k-th wavelet coefficient
		beta = zeros(1,N/(2^j));
		for k=0:(N/(2^j) - 1)	% range of k
			% this the the Psi_{-j,k} th basis element
			Psi = shift(B(j,:),(2^j)*k);
			% dot product with the deconvolution, the estimated coefficient
			% matlab dot product, non-commutative, u.v = sum{conj(u),v}
			% because the norm-square of each Phi or Psi is N, not 1
			% i.e. dot(Psi,Psi) = N
			beta(k+1) = dot(fft(Psi), fdec)/N;

		end
		w = [w beta];
	%%%%%%%%%% computing the variance of the noise
			%%%%%%%%%% at the j-th level
			%%%%%%%%%% Dividing by N since dot(a,b) = dot(fft(a),fft(b))/N	% check on matlab
			%sigmal(j) = sqrt(sigma^2 * dot( (abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2)/N);
			if (length(sigma) == 1)	% that means it is the noise variance
				sigmal(j) = sqrt( dot( (abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);
			else
				%%% For vector valued sigma
				sigmal(j) = sqrt(dot( (abs(fft(sigma)).^2).*(abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2 )/(N.^2));
			end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	end
	% finally, for Phi{-p,k}
	beta = zeros(1,N/(2^p));
	% deconvolution using specific scaling value
	% for vector valued noise
	[fdec, mult] = multi_fdecwien(fsig, fimp, sigma, fori, scaling(p+1));
	%[fdec, mult] = fdecwien(fsig, fimp, sigma_scalar, scaling(p+1));
	for k=0:(N/(2^p) - 1)
		% this the the Psi_{-j,k} th basis element
		Phi = shift(B(p+1,:),(2^p)*k);
		% dot product with the deconvolution, the estimated coefficient
		% matlab dot product, non-commutative, u.v = sum{conj(u),v}
	% because the norm-square of each Phi or Psi is N, not 1
	% i.e. dot(Psi,Psi) = N
		beta(k+1) = dot(fft(Phi), fdec)/N;

		%%%%%%%%%% computing the variance of the noise
		%%%%%%%%%% at the j-th level
		%%%%%%%%%% Dividing by N since dot(a,b) = dot(fft(a),fft(b))/N	% check on matlab
		%sigmal(p+1) = sqrt((sigma^2) .* dot( (abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2)/N);

		
		if (length(sigma) == 1)	% that means it is the noise variance
			sigmal(p+1) = sqrt((sigma^2) .* dot( (abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2)/N);
		else
			%%% For vector valued sigma
			sigmal(p+1) = sqrt(dot( (abs(fft(sigma)).^2).*(abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2 )/(N^2));
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
	w = [w beta];


	% do the thresholding
	% doing the level-dependent thresholding in the same function for now
	thrvec = sigmal.*rho;

	%plotthr(w,p,thrvec)

	[w, ratiounthres, wnoise] = applythres(w, method, p, thrvec);
	%w = keeplarge(w, 2);


	% store the deconvolved signal from the i-th antenna
	z = iwtrans(w,type,p);



%figure 
%plot(z)
%rele(z, testyori)
% plot before windowing
%figure;
%plotsnr(z);
%scaling
%ratiounthres

