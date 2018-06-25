function  plotlevels(z, type,p)
% wnoise is the dropped part of the original signal

w = wtrans(z, type, p);
N = length(w);
M = getbasismat(type,p,N);

% staring index
starting = 1;
len = N/2;
for k=1:p
	ww = zeros(1,N);	% reset it to zero vector every time
	ending = starting+len-1;
	
	ww(starting:ending) = w(starting:ending);
	subplot(p+1,2,k*2-1)
	zz = iwtrans(ww,type,p);
	plot(zz);
	ylabel(num2str(k));
	xlim([1 N]);

	% for fft
	subplot(p+1,2,k*2)
	plotfanita(zz);
	hold on
	plotfanita(M(k,:)*max(abs(fft(zz)))/max(abs(fft(M(k,:)))));
	xlabel('')
	hold off



	% for the next loop
	len = len/2;
	starting = ending + 1;
end

% for the coarsest level Phi_{-p,k}
ww = zeros(1,N);	% reset it to zero vector every time
ww(starting:N) = w(starting:N);
zz = iwtrans(ww,type,p);
subplot(p+1,2,(p+1)*2-1)
plot(zz);
ylabel(strcat('coarse-',num2str(p)));
xlim([1 N]);


% for fft
subplot(p+1,2,(p+1)*2)
plotfanita(zz);
hold on
plotfanita(M(p+1,:)*max(abs(fft(zz)))/max(abs(fft(M(p+1,:)))));
hold off





