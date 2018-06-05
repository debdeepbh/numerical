% thresholding the coefficients of a wavelet transform
% thr is optional
% Correction: hard/soft are the methods
% but value of the threshold are data dependent
function wt = coefthr(w,p,method,thr)
% to get the qth coefficients, do coeff(w, p, q)
N = length(w);

switch method
case 'univ'
	% verify
	sigmacap = 0.5; 	% find this
	lambda = sigmacap * sqrt(2*log2(N)/N);
	thr = zeros(1,p+1) + lambda;
	wt = applythres(w,p,thr);
	
case 'hard'

case 'soft'
case 'maxiset'
	% from Johnstone et al 2004
	% decay parameter of the impulse response
	nu = 1;	% find
	sigma = 0.5;	% should come from the noise variance
	tau_j = 1;	% find
	eta = sqrt(2);	% or eta dependent
	c_n = log2(N)/N;
case 'manual'
	% for p=4
	%thr = [0.8 1 10 10 5];
	wt = applythres(w,p,thr);


end

