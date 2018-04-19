% wavelet compression
function w = compress(z, type, p, K)
temp = keeplarge(wtrans(z, type, p), K);
w = iwtrans(temp, type, p);
