% plot power spectral density in real units using sample rate and 
function plotpanita(z)
fs = 10e9;	% sample rate in GHz
fs = fs/(10^6);	% x axis in MHz
N = length(z);
axsc = fs/N;	% axis scaling factor
x = (1:floor(N/2))*axsc;
plot((1:N)*axsc,10*log10(abs(fft(z))));
% restrict to Nyquist frequency
xlim([1 N*axsc/2])
xlabel('Frequency in MHz')
ylabel('dB')
