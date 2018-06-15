% get the snr of a signal based on
% peak-to-peak/(2*rms)
function snrval = getsnr(z)
snrval = (max(z) - min(z))/(2* sum(abs(z).^2)/sqrt(length(z)));
