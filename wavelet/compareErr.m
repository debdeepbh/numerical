% plots the relative error in pth stage wavelet compression, with respect to K, in q norm
function compareErr(z, maxK, p, q, type1, type2) 
for K=0:maxK-1
	w1(K+1) = norm(z - compress(z,type1, p, K), q)/norm(z, q);
	w2(K+1) = norm(z - compress(z,type2, p, K), q)/norm(z, q);

	fou(K+1)= norm(z - ifft(keeplarge(fft(z),K)))/norm(z, q);
end
x = 1:maxK;
plot(x, w1, x, w2, x, fou);
legend(type1, type2, 'fou');
axis([0 maxK-1 0 1]);
% plot the main graph as well
figure;
plot(0:length(z)-1, z);
axis([0 511]);

