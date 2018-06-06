% incomplete
%%%%%%%%%%%%%%%%%%%
% thresholding the coefficients of a wavelet transform
% thr is optional
% Correction: hard/soft are the rules
% but value of the threshold are data dependent
function wt = coefthr(w, p, rule, method, thr)
% to get the qth coefficients, do coeff(w, p, q)
N = length(w);

switch method
case 'univ'
	% verify
	sigmacap = 0.5; 	% find this
	lambda = sigmacap * sqrt(2*log2(N)/N);
	thr = zeros(1,p+1) + lambda;

case 'forw'	% let's investigate

case 'maxiset'	% incomplete
	% from Johnstone et al 2004
	% decay parameter of the impulse response
	nu = 1;	% find
	sigma = 0.5;	% should come from the noise variance
	tau_j = 1;	% find
	eta = sqrt(2);	% or eta dependent
	c_n = log2(N)/N;
case 'manual'
	% do nothing, ie use the given threshold value
end

% apply the threshold using given rule
wt = applythres(w,p,rule,thr);
