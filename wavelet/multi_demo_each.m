%% comparison of parital, each and averaged(single) on theoretical data
%% check Wiener deconvolution for the 5th single signal
%plot(real(ifft(multi_fdecwien(fft(wax(5,:)), fft(aximp(5,:)), noiseax(5,:), fft(testyori), 1))))

clear all
multi_data

%%% check the avg of the Wiener deconvolution
%for i = 1:15
%	decW(i,:) = real(ifft(multi_fdecwien(fft(wax(i,:)), fft(aximp(i,:)), noiseax(i,:), fft(testyori), 1)));
%end
%avgW = mean(decW);
%%plot(avgW')
%
%% compute the error
%err_wien = rele(avgW, testyori);

type = 'meyer';
p = 5;
rule = 'soft';
rho  = 1;

% forward partial
alpha_par = multi_getscOri_fpar(wax, aximp, testyori, type, p, noiseax)
z_par = multi_fpar(wax, aximp, testyori, type, p,  noiseax, alpha_par, rho, rule);
err_par = rele(z_par, testyori);

sumall = zeros(1,length(wax));
err_each = zeros(1,15);
for i=1:15
% forward each 
	noisenow = std(noiseax(i,:));
	alpha_one = multi_getscOri_fone(wax(i,:), aximp(i,:), testyori, type, p, noisenow);
	z_one = multi_fone(wax(i,:), aximp(i,:), testyori, type, p, std(noiseax(i,:)), alpha_one, rho, rule);
	err_each(i) = rele(z_one, testyori);
	%subplot(5,3,i)
	%plot(z_one); title(['noiseSD: ' num2str(noisenow) ' err: ' num2str(err_each(i))]);
	sumall = sumall + z_one;
end
avgall = sumall/15;
err_all = rele(avgall, testyori);

% forward one
alpha_one = multi_getscOri_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15));
z_one = multi_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15), alpha_one, rho, rule);
err_one = rele(z_one, testyori);

figure
subplot(311)
plot(avgall); title(['avgall ' num2str(err_all)]); xlim([1 length(wax)]);
subplot(312)
plot(z_one); title(['one ' num2str(err_one)]); xlim([1 length(wax)]);
subplot(313)
plot(z_par); title(['par ' num2str(err_par)]); xlim([1 length(wax)]);

figure
noiselevels = std(noiseax');
plot(noiselevels, err_each, 'og');
hold on
plot(noiselevels, zeros(15,1)+err_all); 
plot(noiselevels, zeros(15,1)+err_one); 
plot(noiselevels, zeros(15,1)+err_par); 
%xlim([1 max(noiselevels)]);
legend('individually deconvolved', 'avg of individually deconvolved','deconvolution of avg of observed signal','partially deconvolved then denoised')
xlabel 'noise SD';
ylabel 'relative error';
hold off



