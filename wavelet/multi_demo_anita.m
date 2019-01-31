%% check Wiener deconvolution for the 5th single signal
%plot(real(ifft(multi_fdecwien(fft(wax(5,:)), fft(aximp(5,:)), noiseax(5,:), fft(testyori), 1))))

clear all

[wax, aximp] = prepsig(5);
wax = wax(:,320:320+1023);
for i =1:15
	noiseax(i,:) = randn([1 1024])*std(wax(i,1024-200:1024));
end
testyori = sum(wax)/15/norm(sum(wax)/15);

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
%rho  = 1;
rho = 1;%[2 2 4.5 3 2 10e10];

% forward one
alpha_one = multi_getscOri_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15))
z_one = multi_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15), alpha_one, rho, rule);
err_one = rele(z_one, testyori);

% forward partial
alpha_par = multi_getscOri_fpar(wax, aximp, testyori, type, p, noiseax)
%alpha_par = multi_getoptsc_anita_par(wax, aximp, testyori, type, p, noiseax, 'bisec')
z_par = multi_fpar(wax, aximp, testyori, type, p,  noiseax, alpha_par, rho, rule);
%for i=1:8
%	z_newori = z_par;
%	alpha_par_new = multi_getscOri_fpar(wax, aximp, z_newori, type, p, noiseax)
%	z_par = multi_fpar(wax, aximp, z_newori, type, p,  noiseax, alpha_par, rho, rule);
%	figure 3
%	subplot(8,1,i)
%	plotsnr(z_par)
%end

err_par = rele(z_par, testyori);


subplot(211)
plot(z_par); title(['par ' num2str(err_par)]);
subplot(212)
plot(z_one); title(['one ' num2str(err_one)]);




