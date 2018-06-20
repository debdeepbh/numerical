close all;
p = 1;
type = 'd10';
z = z(1:1024);
x = 1:1024;
z1 = proj(z,type,p);
subplot(211)
plot(z);
title('signal')

subplot(212)
plot(x, z1, x, z-z1);
legend('projection','difference');
figure
subplot(211)
plot(abs(fft(z)));
title('signal')
subplot(212);
plot(x, abs(fft(z1)), x, abs(fft(z-z1)), [1024/(2^p) 1024/(2^p)], [0 max(abs(fft(z)))])

legend('projection','difference');
