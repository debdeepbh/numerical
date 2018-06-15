% estimates the wavelet coefficients after doing deconvolution in Fourier then wavelet domain, based on supplied scaling values
% scaling should be of length p+1
% type is the wavelet filter type, p-th stage
function [w, ratiounthres, thrvec]  = wienforwd(y, K, type, p, sigma, scaling, rho,method)
% sigma is the standard deviation used in wiener deconvolution

N = length(y);

if length(scaling) == 1
	scaling = zeros(1,p+1) + scaling;
end

% level dependent scaling parameter for wavelet thresholding
%rho = 3;

% pad to make them of the same length
paddedimpulse = zeros(1,N);
%starting = floor((length(signal)-length(impulse))/2);
starting=1;
paddedimpulse(starting:starting+length(K)-1) = K;

fsig = fft(y);
fimp = fft(paddedimpulse);



% first get the wavelet basis elements
% j-th row gives the un-shifted wavelet basis
B = getbasismat(type, p, N);


% using equation (22) from Donoho et al 2004 

% matlab allows empty matrix []
w = [];

% varl stores the level dependent leaked noise variance
varl = zeros(1,p+1);

for j=1:p
	% do deconvolution using specific scaling for the level
	[fdec, mult] = fwiendec(fsig, fimp, sigma, scaling(j));

	% prepare beta(k), the estimate for the k-th wavelet coefficient
	beta = zeros(1,N/(2^j));
	for k=0:(N/(2^j) - 1)	% range of k
		% this the the Psi_{-j,k} th basis element
		Psi = shift(B(j,:),(2^j)*k);
		% dot product with the deconvolution, the estimated coefficient
		% matlab dot product, non-commutative, u.v = sum{conj(u),v}
		beta(k+1) = dot(fft(Psi), fdec);

		%%%%%%%%%% computing the variance of the noise
		%%%%%%%%%% at the j-th level
		sigmal(j) = sqrt(sigma^2 * dot( (abs(fft(Psi))./abs(fimp)).^2, abs(mult).^2)/N);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
	w = [w beta];
end
% finally, for Phi{-p,k}
beta = zeros(1,N/(2^p));
% deconvolution using specific scaling value
[fdec, mult] = fwiendec(fsig, fimp, sigma, scaling(p+1));
for k=0:(N/(2^p) - 1)
	% this the the Psi_{-j,k} th basis element
	Phi = shift(B(p+1,:),(2^p)*k);
	% dot product with the deconvolution, the estimated coefficient
	% matlab dot product, non-commutative, u.v = sum{conj(u),v}
	beta(k+1) = dot(fft(Phi), fdec);

	%%%%%%%%%% computing the variance of the noise
	%%%%%%%%%% at the j-th level
	sigmal(p+1) = sqrt(sigma^2 * dot( (abs(fft(Phi))./abs(fimp)).^2, abs(mult).^2)/N);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
w = [w beta];

% because the norm-square of each Phi or Psi is N, not 1
% i.e. dot(Psi,Psi) = N
 w = w./N;


%%%%%%%%%%%%%%%
% plot before thresholding
%plotcoeffs(w,p)
%print('before','-dpng')

% doing the level-dependent thresholding in the same function for now
thrvec = sigmal.*rho;
[w, ratiounthres, wnoise] = applythres(w, method, p, thrvec);
%w = keeplarge(w, 2);

% plot after 
%plotcoeffs(w,p)

