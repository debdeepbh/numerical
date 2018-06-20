% get the isolated component of the wavelet decomposition upto level p
function prj = isolate(z, type, p)

N = length(z);
x = 1:N;
for i=1:p
	Pz = proj(z, type, i); % projection
	Qz = z - Pz;		% residue

	subplot(p+1,2,i*2-1);
	plot(Qz);
	ylabel(num2str(i));
	xlim([1 N]);

	% for fft
	subplot(p+1,2,i*2);
	plot(x, abs(fft(Qz)), [N/(2^(i)) N/(2^(i))], [0 max(abs(fft(Qz)))]);
	xlim([1 N]);
	
	%updating the loop
	z = Pz;
end
subplot(p+1,2,(p+1)*2-1);
plot(Pz);
ylabel('coarsest');
xlim([1 N]);

% for fft
subplot(p+1,2,(p+1)*2);
plot(x, abs(fft(Pz)), [N/(2^(i+1)) N/(2^(i+1))], [0 max(abs(fft(Pz)))]);
xlim([1 N]);

