% wavelet compression
function w = compress(z, type, p, K)
switch type
case 'fou'
	% in case of fourier compression, p is irrelevant, yet necessary
	w = ifft(keeplarge(fft(z),K));
otherwise
	temp = keeplarge(wtrans(z, type, p), K);
	w = iwtrans(temp, type, p);
end
