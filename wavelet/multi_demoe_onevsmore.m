%% compare one (averaged) antenna vs many when both the impulse responses and noise level are the same
% set the parameters
multi_globals

samplesize = 4;
err_par = err_one = err_all = zeros(1,samplesize);
% repeat the experiment 50 times
for count=1:samplesize

	%clear all
	% generate noise every tim
%% 	same  impulse response for each antenna
	multi_data_same

	% forward partial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('par')
	tic
	[alpha_par, ratio_par, sigl_par] = multi_getscOri_fpar(wax, aximp, testyori, type, p, noiseax);
	%alpha_par= multi_getscOri_fpar_plot(wax, aximp, testyori, type, p, noiseax);
	z_par = multi_fpar(wax, aximp, testyori, noiseax, alpha_par);
	toc
	err_par(count) = rele(z_par, testyori);

	% foward individual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('single')
	sumall = zeros(1,length(wax));
	storeall = zeros(15, length(wax));
	for i=1:15
	%i=1;	%change: for one
	% forward each 
		noisenow = std(noiseax(i,:));
		[alpha_single, ratio_single, sigl_single] = multi_getscOri_fone(wax(i,:), aximp(i,:), testyori, type, p, noisenow);
		%alpha_single = multi_getscOri_fone_plot(wax(i,:), aximp(i,:), testyori, type, p, noisenow);
		z_single = multi_fone(wax(i,:), aximp(i,:), testyori, std(noiseax(i,:)), alpha_single);
		%subplot(5,3,i)
		%plot(z_one); title(['noiseSD: ' num2str(noisenow) ' err: ' num2str(err_each(i))]);
		sumall = sumall + z_single;
		storeall(i,:) = z_single;
		inde(i) = rele(z_single, testyori);
	end
	avgall = sumall/15;
	err_all(count) = rele(avgall, testyori);

	% forward one %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('one')
	tic
	[alpha_one, ratio_one, sigl_one] = multi_getscOri_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15));
	%alpha_one = multi_getscOri_fone_plot(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15));
	z_one = multi_fone(sum(wax)/15, sum(aximp)/15, testyori, std(sum(noiseax)/15), alpha_one);
	err_one(count) = rele(z_one, testyori);
	toc

	figure 1;
	plot(z_par);
	xlim([1 length(wax)])
	hold on;
end

xl=1:samplesize;
figure 2;
plot(xl, log(err_all),xl, log(err_one), xl,log(err_par))
legend('avg of individual ForWaRD', 'ForWaRD on avg', 'generalized ForWaRD');
xlim([1 samplesize])
xlabel 'sample'
ylabel 'log of relative error'


%alpha_par
%alpha_single
