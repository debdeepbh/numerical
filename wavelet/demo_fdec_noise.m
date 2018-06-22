% compare errors in different Fourier based deconvolution methods when noise-levels vary
csig = testconv;

smax = 50;

x = [testvec zeros(1,512)];
nx = norm(x);
ea =zeros(1,smax);
ew=zeros(1,smax);

snra = snrw = zeros(1,smax);

for sigma=1:smax
	ns = randn([1 1024])*sigma; 
	sig = csig + ns;

	% error in allpass
	siga = decall(sig, K, 'allp', 1);
	sigw = decwien(sig, K, randn([1 1024])*sigma, 1);
	ea(sigma) = norm(siga - x)/nx;
	ew(sigma) = norm(sigw - x)/nx;

	snra(sigma) = getsnr(siga);
	snrw(sigma) = getsnr(sigw);

	if (sigma == 20)
		figure(3);
		plotsnr(sigw);
	end
end
figure(1);
plot(1:smax,log10(ea),1:smax,log10(ew));
legend('allpass','Wiener');
title('relative error in 2-norm');
xlim([1 smax]);
xlabel('noise standard deviation')
ylabel('log10 of relative l-2 error')

figure(2);
plot(1:smax,log10(snra),1:smax,log10(snrw));
legend('allpass','Wiener');
title('SNR using peak-to-peak/(2*rms(outside))');
xlim([1 smax]);
xlabel('noise standard deviation')
ylabel('log10(snr)')



