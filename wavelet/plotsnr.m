function plotsnr(z)
[snrval, M, m, a, b] = getsnr(z);
	plot(z)
	hold on;
	plot([a a],[m M]);
	plot([b b],[m M]);
	xlim([1 length(z)]);
	title(num2str(snrval));
	hold off;

