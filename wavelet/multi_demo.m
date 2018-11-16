%% check Wiener deconvolution for the 5th single signal
%plot(real(ifft(multi_fdecwien(fft(wax(5,:)), fft(aximp(5,:)), noiseax(5,:), fft(testyori), 1))))

multi_data

%% check the avg of the Wiener deconvolution
for i = 1:15
	decW(i,:) = real(ifft(multi_fdecwien(fft(wax(i,:)), fft(aximp(i,:)), noiseax(i,:), fft(testyori), 1)));
end
avgW = mean(decW);
plot(avgW')

% compute the error
disp 'Wiener'
rele(avgW, testyori)

% forward partial
alphanew = multi_getscOri_fpar(wax, aximp, testyori, 'meyer', 5, noiseax);
z = multi_fpar(wax, aximp, testyori, 'meyer', 5, noiseax, alphanew, 1, 'soft');
disp 'fpar'
rele(z, testyori)

% forward one
betanew = multi_getscOri_fone(sum(wax)/15, sum(aximp)/15, testyori, 'meyer', 5, std(sum(noiseax)/15));
z1 = multi_fone(sum(wax)/15, sum(aximp)/15, testyori, 'meyer', 5, std(sum(noiseax)/15), betanew, 1, 'soft');
disp 'fone'
rele(z1, testyori)




