% generate a threshold vector of length p+1 for a p-th stage
% wavelet transform
% supply the variance, manual threshold optional
function thr = genthr(p, var, method, manual)

switch method
case 'univ'	% universal thresholding using manual value
	lambda = var * sqrt(2*log2(N)/N);
	thr = zeros(1,p+1) + lambda;
case 'maxis'	% maxiset thresholding

case 'man'	% supply manually, a vector
	thr = manual;
end


