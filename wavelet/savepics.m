% save the plots from each run
clear all;
getdata; 	%just to get the noise data, although would like to have the noise data during the time of recording
for antindex=1:20
	[ax, aximp] = prepsig(antindex);
	subplot(5,4,antindex);
	wienforwdant;
	%finalfile = strcat('sig_',num2str(antindex));
	%demo_save(1,finalfile, '-S800,200');
end

