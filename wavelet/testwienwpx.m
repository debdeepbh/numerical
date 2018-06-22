% compare SNR for Fourier based deconvolutions
for i=1:8
	z = wpx(i,:)(1:1024);

	figure(1);
	subplot(4,2,i);
	plotsnr(z);

	figure(2);
	subplot(4,2,i);
	plotsnr(decwien(z,K,noise,1));

	figure(3);
	subplot(4,2,i);
	plotsnr(decall(z,K,'allp',1));
end

