i=1;
rhoarr = 1:0.5:10;
snv = zeros(1,length(rhoarr));
err = zeros(1,length(rhoarr));

x = [testvec zeros(1,512)];	%original signal

for rho=rhoarr
	z = testy;
	type = 'd10';
	p = 5;
	sigma =5;
	alpha = sc;
	%alpha = [0.325 0.275 0.375 0.425 0.725 0.675];
	method = 'soft';
[w, ratiounthr, thrvec] = wienforwd(z(1:1024),K(1:1024),type,p,sigma,alpha, rho,method);
 z = iwtrans(w, 'd10', p);
 err(i) = norm(x-z)/norm(x);
 snv(i) = getsnr(z);
 i = i+1;
% subplot(4,1,rho)
% plot(z)
 end
figure(1)
plot(rhoarr, snv)
% figure(2)
% plot(rhoarr, err)

%Hard & 0.10177 & 0.10818 & 0.10995 & 0.11546
%Soft & 0.10518 & 0.10882 & 0.11030 & 0.11532 
%Wiener & 0.086166 & 0.094167 & 0.096968 & 0.10215 
