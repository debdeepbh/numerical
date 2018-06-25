function plotsnr(z)
[snrval, M, m, a, b] = getsnr(z);
	plot(z)
	hold on;
	plot([a a],[m M],'r');
	plot([b b],[m M],'r');
	xlim([1 length(z)]);
	title(num2str(snrval));
	hold off;

