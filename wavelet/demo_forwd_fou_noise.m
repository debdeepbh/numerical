% compare errors between Wiener and Forward method
csig = testconv;

smax = 40;
avgc = 5;	% avg over 5 random draws

type = 'meyer';
p= 5;
rho = [1 1 2 1 1 10e10];	% change this if does not work
%rho = 1;
method = 'soft';

sigmalist = 5:5:smax;

sl = length(sigmalist);

ef = ew = ea = zeros(avgc,sl);

snrf = snrw =  zeros(avgc,sl);
sc_temp = zeros(avgc,sl,p+1);

i=1;
for sigma=sigmalist
%sigma=5;

	%for j=1:avgc
	j=1;
		ns = randn([1 1024])*sigma; 
		sig = testconv + ns;
		%sc_temp(j,i,:) = getoptsc(sig,K,type,p,ns,'bisec');

		sigf = iwtrans(wienforwd(sig, K, type, p, ns, sc_temp(j,i,:), rho, method),type,p);
		sigf = proj(sigf,type,2);

		sigw = decwien(sig, K,ns, 1);
		siga = decall(sig, K,'allp', 1);

		ef(j,i) = rele(sigf,testyori);
		ew(j,i) = rele(sigw,testyori);
		ea(j,i) = rele(siga,testyori);

		snrf(j,i) = getsnr(sigf);
		snrw(j,i) = getsnr(sigw);
		snra(j,i) = getsnr(siga);
	%end
	%aef(sigma) = median(ef(:,sigma));
	%aew(sigma) = median(ew(:,sigma));
	%asnrf(sigma) = median(snrf(:,sigma));
	%asnrw(sigma) = median(snrw(:,sigma));


	%%if (sigma == 20)
	%	figure(3);
	%	subplot(sl,2,i*2-1)
	%	plotsnr(sigf);
	%	subplot(sl,2,i*2)
	%	plotsnr(sigw);

	%%end
	i=i+1;
end

length(ef(1,:))
length(sigmalist)
figure(1);
%plot(1:smax,log10(aef),1:smax,log10(aew));
plot(sigmalist,log10(ef(j,:)),sigmalist,log10(ew(j,:)),sigmalist,log10(ea(j,:)));
legend('forward','Wiener','allpass');
title('log10 of relative error in 2-norm');
xlim([5 smax]);
xlabel('noise standard deviation')
ylabel('relative l-2 error')

figure(2);
%plot(1:smax,log10(asnrf),1:smax,log10(asnrw));
plot(sigmalist,snrf(j,:),sigmalist,snrw(j,:),sigmalist,snra(j,:));
legend('forward','Wiener','allpass');
title('SNR using peak-to-peak/(2*rms(outside))');
xlim([5 smax]);
xlabel('noise standard deviation')
ylabel('snr')



