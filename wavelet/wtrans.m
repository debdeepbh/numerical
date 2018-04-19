%wavelet transform, general
function w = wtrans(z, type, p)
N = length(z);
[u, v] = filt(type, N);

util = ifft(conj(fft(u)));
vtil = ifft(conj(fft(v)));
 
w = wrec(z, N/2^(p-1), util, vtil);

