% plots the relative error in pth stage wavelet compression, with respect to K, in q norm
function relerr(z, maxK, p, q, ef, tf) 
for K=0:maxK-1
	w1(K+1) = norm(z - compress(z, 'd2', p, K), q)/norm(z, q);
	w2(K+1) = norm(z - compress(z, 'd4', p, K), q)/norm(z, q);
	w3(K+1) = norm(z - compress(z, 'd6', p, K), q)/norm(z, q);
	w4(K+1) = norm(z - compress(z, 'd10', p, K), q)/norm(z, q);

	fou(K+1)= norm(z - ifft(keeplarge(fft(z),K)))/norm(z, q);
end
x = 1:maxK;
plot(x, w1, x, w2, x, w3, x, w4, x, fou);
legend('d2','d4','d6','d10', 'fou');
axis([0 maxK-1 0 1]);
print(ef,'-dpng');
% plot the main graph as well
figure;
plot(0:length(z)-1, z);
axis([0 511]);
print(tf,'-dpng');

