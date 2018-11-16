% streamlining the deconvolution while deconvoling each antenna. Wavelet thresholding done in the end
% estimates the wavelet coefficients after doing deconvolution in Fourier then wavelet domain, based on supplied scaling values
% scaling should be of length p+1
% type is the wavelet filter type, p-th stage
%function [w, ratiounthres, thrvec]  = wienforwd(y, K, type, p, sigma, scaling, rho,method)
% sigma is the standard deviation used in wiener deconvolution


type = 'meyer';
p = 5;
method = 'soft';
rho = 1;

noisearr = std(noiseax(:,1:1024)')';


zaxf = zaxa = zaxw = zeros(size(wax));

snrval = zeros(1,15);

% compute each time or once
%sc_ax = zeros(15,p+1)+1;
sc_ax = load('sc_ax');


N = 1024;
M = 15;

% get the wavelet basis
B = getbasismat(type, p, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet Shrinkage

sigmal = zeros(M,p+1);
sigmalavg = zeros(1,p+1);

% using equation (22) from Donoho et al 2004 

% all the antennas
for i =1:M
	fsig = fft(wax(i,:));
	fimp = fft(aximp(i,:));
	fori = fft(testyori);
	sigma = noisearr(i);	% using scalar noise, i.e standard deviation

	% to be chosen later
	scaling = sum(sc_ax)/M;

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
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end
		w = [w beta];
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
			sigmal(i,p+1) = sqrt((sigma^2) .* dot( (abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2)/N);
		else
			%%% For vector valued sigma
			sigmal(i,p+1) = sqrt(dot( (abs(fft(sigma)).^2).*(abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2 )/(N^2));
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
	w = [w beta];

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

%%%%%%%%% windowing to cut off higher frequencies
%z = windowfreq(z, 0, 180/5000, 600/5000, 0.75);
%
%%figure;
%plotsnr(z);



%figure;
%plotpanita(z);

%isolate(z,type,p);

% plot after 
%plotcoeffs(w,p)

scaling
ratiounthres

