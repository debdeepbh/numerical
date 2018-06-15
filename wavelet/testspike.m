% test the values of forward
function testspike(z, type, p, sigma,  alpha, rho, method)
% example
% plot(iwtrans(decfw(wex(3,:)(1:1024),K(1:1024),'d6',4,'wien',[10 400 50 50 50], 4),'d6',4))
% testfwd(3,'d6',4,[10 400 50 50 50],10)
% optimal match at testfwd(3,'d6',4,[110 100 144 432 320],2)
% Strategy: first get the optimal values using rho = 1
% then adjust rho looking at the wavelet transform
% testspike(testy,'d6',4,5,[0.35 0.35 0.3 0.65 0.7],[2 2 1.5 1 1],'soft')
% For anita data: testspike(wex(3,:), 'd10', 5,sqrt(8.5),[0.3 0.3 0.36 0.64 0.49 4],[4 4 4 1 1 1],'hard')
% bumpless noiseless barebone: testspike(wex(5,:), 'd10', 5,sqrt(8.5),[0.3 0.3 0.36 0.64 0.49 4],[4 4 4 1 3 3],'hard')

getdata;
close all;
[w, ratiounthr, thrvec] = wienforwd(z(1:1024),K(1:1024),type,p,sigma,alpha, rho,method);

plotthr(w,p,thrvec);
%plotcoeffs(w,p);

figure;
subplot(211)
plot(z)
xlim([1 length(w)])
switch type
case 'bior13d'
	type = 'bior13d'
end
subplot(212)
plot(iwtrans(w,type,p));
xlim([1 length(w)])


figure;
plot(log10(abs(fft(iwtrans(w,type,p))).^2));
title('Power Spectrum of output in log10 scale');
xlim([1 length(w)])

% print the threshold proportion
ratiounthr
