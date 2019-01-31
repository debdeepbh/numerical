% plot the mse with all values of parameters
tic

paramvals = 0.001:0.1:1;
errVals = zeros(size(paramvals));
alpha = zeros(1, p+1);
bestparam = zeros(1, p+1);
minerr = zeros(1, p+1);
for j=1:p+1
	i = 1;
	for param=paramvals
		alpha(j) = param;
		z_par = multi_fpar(wax, aximp, testyori, noiseax, alpha);
		z_one = multi_fone(mean(wax), mean(aximp), testyori,
		std(mean(noiseax)), alpha);
		errVals(i) = rele(z_par, testyori);
		one_errVals(i) = rele(z_one, testyori);
		i = i+1;
	end
	[minerr(j) paramInd] = min(errVals);
	[one_minerr(j) one_paramInd] = min(one_errVals);
	bestparam = paramvals(paramInd(1));
	one_bestparam = paramvals(one_paramInd(1));
	alpha(j) = bestparam;
	one_alpha(j) = one_bestparam;
	fprintf('%d\n',j);
	fflush(stdout)
end

toc


%%%  Print things
alpha
minerr
one_alpha
one_minerr
subplot(211)
		z_par = multi_fpar(wax, aximp, testyori,  noiseax, alpha);
		plot(z_par)
		z_one = multi_fone(mean(wax), mean(aximp), testyori,
		 std(mean(noiseax)), one_alpha);
subplot(212)
		plot(z_one)


