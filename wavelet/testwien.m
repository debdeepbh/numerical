% test wiener deconvolution
clear all;
getdata;

alpha = 1;	%the scaling parameter


sig = testy(1:1024);
imp = K(1:1024);
fsig = fft(sig);
fimp = fft(imp);

naive = fsig./fimp;

sigma = 5;
hsq = abs(fimp).^2;

fori = abs(fft([testvec zeros(1,512)])).^2;	% hiller-chin justifies this choice
%fori = hsq;

mult = (hsq./(hsq + alpha*1024*sigma^2./fori));
%mult2 = (hsq./(hsq + alpha*1024*(abs())./fori));

subplot(211)
plot(testvec);
subplot(212)
plot(real(ifft(mult.*naive)))
