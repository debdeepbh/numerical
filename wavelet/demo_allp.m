h = gcf;
set(h,'PaperPositionMode','auto');  

subplot(221)
plot(testy);
xlim([0 1024]);
title('Observed noisy signal');

z = wex(3,:);
subplot(222)
plot(z);
xlim([0 1024]);
title('Observed noisy signal');

subplot(223)
plot(decall(testy,K,'allp',1));
xlim([0 1024]);
title('Deconvolved with allpass');

subplot(224)
plot(decall(z,K,'allp',1));
xlim([0 1024]);
title('Deconvolved with allpass');

print gcf -dpdf '-S800,500' '/home/debdeep/gdrive/anita/frankdata/pic/demo_allp.pdf'
