% plot the signals used in the latex file
h = gcf;
set(h,'PaperPositionMode','auto');  

i=1;
subplot(2,2,i); i = i+1;
plot([testvec zeros(1,512)]);
xlim([1 1024]);
title('Orignal signal: x');

subplot(2,2,i); i = i+1;
plot(K);
xlim([1 1024]);
title('Impulse response: h');

xconvo = conv(testvec,K);
subplot(2,2,i); i = i+1;
plot(xconvo);
xlim([1 1024]);
title('Convolved signal: x*h')


subplot(2,2,i); i = i+1;
plot(testy);
xlim([1 1024]);
title('Observed (noisy) signal: y = x*h+n');

%set(ha(1:p),'XTickLabel',''); 

%set(ha,'YTickLabel','');

print gcf -dpdf '-S800,400' '/home/debdeep/gdrive/anita/frankdata/pic/demo_signal.pdf'

