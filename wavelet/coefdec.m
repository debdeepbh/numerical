% compute the wavelet coefficients based on observed signal y,
% impulse response K and p-th stage wavelet decomposition using
% a type of wavelets
% assumption: y, K are dyadic of the same length
function w = coefdec(y, K, type, p)
N = length(y);
% first get the wavelet basis elements
% j-th row gives the un-shifted wavelet basis
B = getbasismat(type, p, N);

% deconvolution method here, on the fourier side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dec = fft(y)./fft(K); 	% naive
dec = fft(y)./fft(K).*abs(fft(K)); % allpass
%%%%%%% other methods should go here %%%%%%%

% using equation (22) from Donoho et al 2004 

% matalb allows empty matrix []
w = [];

for j=1:p
	% prepare beta(k), the estimate for the k-th wavelet coefficient
	beta = zeros(1,N/(2^j));
	for k=0:(N/(2^j) - 1)	% range of k
		% this the the Psi_{-j,k} th basis element
		Psi = shift(B(j,:),(2^j)*k);
		% dot product with the deconvolution, the estimated coefficient
		% matlab dot product, non-commutative, u.v = sum{conj(u),v}
		beta(k+1) = dot(fft(Psi), dec);
	end
	w = [w beta];
end
% finally, for Phi{-p,k}
beta = zeros(1,N/(2^p));
for k=0:(N/(2^p) - 1)
	% this the the Psi_{-j,k} th basis element
	Phi = shift(B(p+1,:),(2^p)*k);
	% dot product with the deconvolution, the estimated coefficient
	% matlab dot product, non-commutative, u.v = sum{conj(u),v}
	beta(k+1) = dot(fft(Phi), dec);
end
w = [w beta];

% for some reason, (probably the normalization while applying plancheral)
 w = w./N;


