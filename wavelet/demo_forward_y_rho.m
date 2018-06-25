% plot the result of threhsolding
rholist = 1:5;
eh = es = zeros(1,length(rholist));
i=1;
for rho=rholist
	subplot(2,5,i);
	zh = iwtrans(wienforwd(testy, K, 'meyer', 5, testnoise, sc_y_m, rho, 'hard'),'meyer',5);
	plot(zh);
	title(strcat('rho=',num2str(rho)));
	xlim([1 1024])
	i=i+1;
end
for rho=rholist
	subplot(2,5,i);
	zh = iwtrans(wienforwd(testy, K, 'meyer', 5, testnoise, sc_y_m, rho, 'soft'),'meyer',5);
	plot(zh);
	title(strcat('rho=',num2str(rho)));
	xlim([1 1024])
	i=i+1;
end
