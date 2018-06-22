% get threshold values based on noise data, after Wiener deconvolution
function getthrfromnoise(noise, K, type, p)
% strange because this is where you use the information from noise as well
out = decall(noise, K, 'wien', 1);
w = wtrans(out, type, p);
for i=1:p+1
	coefvec = coeff(w, p, i);
	% picking the maximum value as the threshold
	thr(i) = max(abs(coefvec));
end

