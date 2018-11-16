%% check Wiener deconvolution for the 5th single signal
%plot(real(ifft(multi_fdecwien(fft(wax(5,:)), fft(aximp(5,:)), noiseax(5,:), fft(testyori), 1))))

clear all
multi_data

%% check the avg of the Wiener deconvolution
for i = 1:15
	decW(i,:) = real(ifft(multi_fdecwien(fft(wax(i,:)), fft(aximp(i,:)), noiseax(i,:), fft(testyori), 1)));
end
avgW = mean(decW);
%plot(avgW')

% compute the error
err_wien = rele(avgW, testyori);

type = 'd10';
p = 5;
rule = 'soft';
rho  = 1;

% forward partial
alpha_par = multi_getscOri_fpar(wax, aximp, testyori, type, p, noiseax)
z_par = multi_fpar(wax, aximp, testyori, type, p,  noiseax, alpha_par, rho, rule);
err_par = rele(z_par, testyori);

% forward one
alpha_one = multi_getscOri_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15))
z_one = multi_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15), alpha_one, rho, rule);
err_one = rele(z_one, testyori);

subplot(211)
plot(z_par); title(['par ' num2str(err_par)]);
subplot(212)
plot(z_one); title(['one ' num2str(err_one)]);




