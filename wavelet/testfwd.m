% test the values of forward
function ratiothr = testfwd(anindex, type, p, alpha, rho)
% example
% plot(iwtrans(decfw(wex(3,:)(1:1024),K(1:1024),'d6',4,'wien',[10 400 50 50 50], 4),'d6',4))
% testfwd(3,'d6',4,[10 400 50 50 50],10)
% optimal match at testfwd(3,'d6',4,[110 100 144 432 320],2)
getdata;
close all;

figure;
subplot(211)
plot(wex(anindex,:)(1:1024))

subplot(212)
[w, ratiothr] = decfw(wex(anindex,:)(1:1024),K(1:1024),type,p,'wien',alpha, rho);
plot(iwtrans(w,type,p));
ratiothr = ratiothr*1024;
