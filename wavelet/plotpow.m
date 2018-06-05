% plot the power of a signal in fourier domain
function plotpow(sig)
plot(log10(abs(fft(sig)).^2));
% get the name of the variable as string
name = inputname(1);
xlim([1 length(sig)/2]);
titlestr = strcat('log_{10}(|fft(',name,')|^2)');
title(titlestr);

