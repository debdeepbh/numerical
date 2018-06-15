% different deconvolution methods
function w = decall(signal, impulse, method, scaling)
% scaling is optional
% var is also tau for tikhonov

% figure out optimal padding
% pad impulse by zero to make them of the same length
% assuming here that impulse is of smaller length than signal
N = length(signal);
paddedimpulse = zeros(1,N);
%starting = floor((length(signal)-length(impulse))/2);
starting=1;
paddedimpulse(starting:starting+length(impulse)-1) = impulse;

% do the deconvolution on the Fourier side, and then ifft
w = real(ifft(fdecall(fft(signal), fft(paddedimpulse), method, scaling)));
    


