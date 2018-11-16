function [z, ratiounthres] = multi_fpar(wax, aximp, testyori, type, p, noiseax, sc_ax, rho, method)

noisearr = std(noiseax')';


 [M, N] = size(wax);
% get the wavelet basis
B = getbasismat(type, p, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet Shrinkage

sigmal = zeros(M,p+1);
sigmalavg = zeros(1,p+1);

% using equation (22) from Donoho et al 2004 

fori = fft(testyori);
% all the antennas
for i =1:M
	fsig = fft(wax(i,:));
	fimp = fft(aximp(i,:));
	sigma = noisearr(i);	% using scalar noise, i.e standard deviation

	% to be chosen later
	scaling = sc_ax;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% estimate the wavelet coefficients 
	%% and the noise variance of summed data

	% matlab allows empty matrix []
	w = [];

	% varl stores the level dependent leaked noise variance
	varl = zeros(1,p+1);

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
					sigmal(i,j) = sqrt( dot( (abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2 * sigma^2 )/N);
				else
					%%% For vector valued sigma
					sigmal(i,j) = sqrt(dot( (abs(fft(sigma)).^2).*(abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2 )/(N.^2));
				end

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

	end
	w = [w beta];
	%%%%%%%%%% computing the variance of the noise %%%%%%%%%% at the j-th level
	%%%%%%%%%% Dividing by N since dot(a,b) = dot(fft(a),fft(b))/N	% check on matlab
	%sigmal(p+1) = sqrt((sigma^2) .* dot( (abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2)/N);

	
	if (length(sigma) == 1)	% that means it is the noise variance
		sigmal(i,p+1) = sqrt((sigma^2) .* dot( (abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2)/N);
	else
		%%% For vector valued sigma
		sigmal(i,p+1) = sqrt(dot( (abs(fft(sigma)).^2).*(abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2 )/(N^2));
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% store the scaled-wiener deconvolved signal's wavelet transform for i-th antenna
	antWienw(i,:) = w;

end

% get the  avg of the antennas, wtrans 
avgWienw = sum(antWienw)/M;

% get the noise sd for avg of antenna
for j=1:p+1
	%sigmalavg(j) = sqrt(sum((sigmal(:,j)).^2)/M);	% avg over i %probably wrong M
	sigmalavg(j) = sqrt(sum((sigmal(:,j)).^2))/M;	% avg over i
end


%%%%%%%%%%%%%%%
% plot before thresholding
%plotcoeffs(w,p)
%print('before','-dpng')

% doing the level-dependent thresholding in the same function for now
thrvec = sigmalavg.*rho;

%%% print before applying threshold
%figure;
%plot(iwtrans(w,type,p))
%%plotthr(w,p,thrvec);
%title('before applying threshold')

plotthr(avgWienw, p, thrvec)


[w, ratiounthres, wnoise] = applythres(avgWienw, method, p, thrvec);
%w = keeplarge(w, 2);

% print the threshold values
%thrvec

% plot the threshold values after thresholding
%plotthr(w,p,thrvec);

z = iwtrans(w,type,p);

figure 2
plot(z)
rele(z, testyori)
% plot before windowing
%figure;
%plotsnr(z);

scaling
ratiounthres
