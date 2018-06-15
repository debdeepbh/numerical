% iterative
% test wiener deconvolution
clear all;
getdata;

maxiter = 7; 
alpha = 1;	%the scaling parameter

sig = testy(1:1024);
imp = K(1:1024);
fsig = fft(sig);
fimp = fft(imp);

naive = fsig./fimp;

sigma = 5;
hsq = abs(fimp).^2;

% first choose the observed signal as a crude choice of desired signal
fori = abs(fsig);	% hiller-chin justifies this choice
%fori = hsq;
fnoisesq = sigma^2;

figure(1);
for j=1:maxiter
	mult = (hsq./(hsq + alpha*1024*fnoisesq./(fori.^2)));
	fout = mult.*naive;
	
	% plotting
	z = real(ifft(fout));
	%if j==maxiter
		subplot(maxiter,1,j);
		plot(z);
		title(num2str(j));
	%end

	% update
	fori = fout;
	% correction term from Hiller-Chin for convergence
	fcorrec = (sigma^2).*(fori.^2)./(hsq.*(fori.^2) + sigma^2);
	fori = fori + fcorrec;
	%fnoisesq = fsig.^2 - (fori.^2).*hsq;
	%fnoisesq = fnoisesq*0.8;
end

