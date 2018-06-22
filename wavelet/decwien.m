% different deconvolution methods
function w = decwien(signal, impulse, sigma, scaling)
% scaling is optional
% sigma is either raw noise of noise variance

% figure out optimal padding
% pad impulse by zero to make them of the same length
% assuming here that impulse is of smaller length than signal
N = length(signal);
paddedimpulse = zeros(1,N);
%starting = floor((length(signal)-length(impulse))/2);
starting=1;
paddedimpulse(starting:starting+length(impulse)-1) = impulse;

% do the deconvolution on the Fourier side, and then ifft
w = real(ifft(fdecwien(fft(signal), fft(paddedimpulse), sigma, scaling)));
    


