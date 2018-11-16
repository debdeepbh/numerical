%% check Wiener deconvolution for the 5th single signal
%plot(real(ifft(multi_fdecwien(fft(wax(5,:)), fft(aximp(5,:)), noiseax(5,:), fft(testyori), 1))))

%% check the avg of the Wiener deconvolution
for i = 1:15
	decW(i,:) = real(ifft(multi_fdecwien(fft(wax(i,:)), fft(aximp(i,:)), noiseax(i,:), fft(testyori), 1)));
end
avgW = mean(decW);
plot(avgW')


% compute the error
rele(avgW, testyori)
