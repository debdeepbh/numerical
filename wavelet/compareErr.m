% plots the relative error in pth stage wavelet compression, with respect to K, in q norm
function compareErr(z, typelist, maxK, p, normlist)
% e.g. typelist = {'shan'; 'db5'}
% e.g. normlist = {1,2,inf}

% number of norms to compare
n = length(normlist);

x = 1:maxK;
% plot the data itself, spanning the first row
subplot(2,n,[1:n]);
plot(0:length(z)-1, z);
axis([0 length(z)]);
title('signal')

% 1 plot each for each norm
for nor=1:n
	q = normlist{nor};
	x = 1:maxK;
	subplot(2,n,n+nor);
	plot(x,relerr(z,typelist,maxK,p,q))
	axis([0 maxK-1 0 1]);
	legend(typelist);
	title([ num2str(q), '-norm']);
	ylabel(['relative error']);
	xlabel(['size of compressed signal']);
end


%print(tf,'-dpng');
