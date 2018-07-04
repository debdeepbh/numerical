% deconvolve the antennas then add them up
type = 'meyer';
p = 5;
% noise is loaded from another file
%noiseax = nx(:,320:(320+1023));
noiseax = zeros(15,1)+11.2; %nx(:,1:(1+1023));
%noiseax = 18.2;
%rho = 2;
%rho =2; \
rho = [10e10 10e10 4 3 2 10e10];

%wax = ax(:,320:(320+1023));
wax = ax(:,420:(420+1023));

zaxf = zaxa = zaxw = zeros(size(wax));

snrval = zeros(1,15);

% compute each time or once
%sc_ax = zeros(15,p+1);
sc_ax = load('sc_ax');

for i=1:15
	K = aximp(i,:);
	%if (i==1)
	%	sc_ax(i,:) = getoptsc(wax(i,:),K,type,p,noiseax(i,:), 'bisec');
	%end
	ww = wienforwd(wax(i,:),K, type, p, noiseax(i,:), sc_ax(i,:), rho, 'soft');
	zaxf(i,:) = iwtrans(ww,type,p);
	snrval(i) = getsnr(zaxf(i,:));

	% trying this
	% shannon, keeplarge 3 works
	%zaxf(i,:) = iwtrans(keeplarge(ww,3),type,p);

	zaxa(i,:) = decall(wax(i,:),K,'allp',1);
	zaxw(i,:) = decwien(wax(i,:),K,noiseax(i,:),1);

end


% delete all
close all

figure;
plot(zaxf');
xlim([0 1024])

figure;
%sumz = sum(zaxf.*(snrval'));
sumz = sum(zaxf)/15;
%sumz = median(zaxf);	% median stacking
sumzdn = sumz;
%sumzdn = compress(sumz,'meyer',3,20);

plotsnr(sumzdn);
ylabel('forwd then sum');

figure;
plotsnr(sum(zaxa)/15);
ylabel('allp then sum');

figure
plotsnr(decall(sum(wax)/15,sum(aximp)/15,'allp',1));
ylabel('sum then allp')


figure;
plotsnr(sum(zaxw)/15);
ylabel('wiener then sum');

figure
plotsnr(decwien(sum(wax)/15,sum(aximp)/15,sum(noiseax)/15,1));
ylabel('sum then wiener')




	


