%% re-generate the noise many times and plot the relative error on theoretical data 
% set the parameters
type = 'meyer';
p = 5;
rule = 'soft';
rho  = 1;

samplesize = 50;
err_par = err_one = err_all = zeros(1,samplesize);
% repeat the experiment 50 times
for count=1:samplesize

	%clear all
	% generate noise every time
	multi_data
	% forward partial
	alpha_par = multi_getscOri_fpar(wax, aximp, testyori, type, p, noiseax);
	z_par = multi_fpar(wax, aximp, testyori, type, p,  noiseax, alpha_par, rho, rule);
	err_par(count) = rele(z_par, testyori);

	sumall = zeros(1,length(wax));
	for i=1:15
	% forward each 
		noisenow = std(noiseax(i,:));
		alpha_one = multi_getscOri_fone(wax(i,:), aximp(i,:), testyori, type, p, noisenow);
		z_one = multi_fone(wax(i,:), aximp(i,:), testyori, type, p, std(noiseax(i,:)), alpha_one, rho, rule);
		%subplot(5,3,i)
		%plot(z_one); title(['noiseSD: ' num2str(noisenow) ' err: ' num2str(err_each(i))]);
		sumall = sumall + z_one;
	end
	avgall = sumall/15;
	err_all(count) = rele(avgall, testyori);

	% forward one
	alpha_one = multi_getscOri_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15));
	z_one = multi_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15), alpha_one, rho, rule);
	err_one(count) = rele(z_one, testyori);

end


xl=1:samplesize;
plot(xl, err_all,xl, err_one, xl,err_par)
legend('avg of induvidually deconvolved', 'deconvolved avg', 'partially dec wavelet denoise');
xlim 'sample'
ylim 'relative error'



