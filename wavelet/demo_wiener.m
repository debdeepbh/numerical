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
plot(decall(testy,K,'wien',1));
xlim([0 1024]);
title('Deconvolved with Wiener');

subplot(224)
plot(decall(z,K,'wien',1));
xlim([0 1024]);
title('Deconvolved with Wiener');

print gcf -dpdf '-S800,500' '/home/debdeep/gdrive/anita/frankdata/pic/demo_wien.pdf'
