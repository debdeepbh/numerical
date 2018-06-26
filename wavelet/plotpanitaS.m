% plot power spectral density in real units using sample rate and 
function plotpanitaS(z)
fs = 2.6e9;	% sample rate in GHz
fs = fs/(10^6);	% x axis in MHz
N = length(z);
axsc = fs/N;	% axis scaling factor


F = 10*log10(abs(fft(z)));
FF = F(1:N/2+1);	% this is fhat(0) to fhat(N/2)
			% this should be scaled to [0 fs/2]
xaxis = (0:N/2)*axsc;	% this is the axis where it should be plotted
plot(xaxis, FF);
xlim([0 N*axsc/2]);
xlabel('Frequency in MHz')
ylabel('dB')
