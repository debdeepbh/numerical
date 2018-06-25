% compute relative error when rho changes

rholist = 1:0.5:4;
eh = es = zeros(1,length(rholist));
i=1;
for rho=rholist
	zh = iwtrans(wienforwd(testy, K, 'meyer', 5, testnoise, sc_y_m, rho, 'hard'),'meyer',5);
	eh(i) = rele(zh,testyori);
%	eh(i) = getsnr(zh);
	zs = iwtrans(wienforwd(testy, K, 'meyer', 5, testnoise, sc_y_m, rho, 'soft'),'meyer',5);
	es(i) = rele(zs,testyori);
%	es(i) = getsnr(zs);

	i=i+1;
end

figure(1)
plot(rholist,eh)
hold on
plot(rholist,es)
xlabel('rho');
%ylabel('snr');
ylabel('relative error in 2-norm');
%title('SNR for level-independent thresholding parameter')
title('Relative error vs level-independent thresholding parameter')
legend('hard','soft')

