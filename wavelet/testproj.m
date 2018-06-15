testvec_4
for i=0:5
	figure(1);
	subplot(5,1,i+1);
	p = proj(testvec,'d10',i);
	plot(p);
	figure(2);
	subplot(5,1,i+1);
	plot(abs(fft(p)));
end


