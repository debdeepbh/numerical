alpha = zeros(1,p+1) + 1;
for i =1:p+1
	count =1;
	for a = 1e-3:1e-1:1

		alphaval(count) = a;
		alpha(i) = a;
		errvalone(count) = rele(multi_fone(sum(wax)/15, sum(aximp)/15, testyori, type, p, std(sum(noiseax)/15), alpha, rho, rule), testyori);
		errvalpar(count) = rele(multi_fpar(wax,aximp, testyori, type, p, noiseax, alpha, rho, rule), testyori);
		count = count + 1;
	end
	figure 3
	subplot(p+1, 1, i);
	plot(alphaval, log(errvalone))
	figure 4
	subplot(p+1, 1, i);
	plot(alphaval, log(errvalpar))

end
		
