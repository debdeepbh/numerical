% test the values of forward
function getthr(anindex, type, p, alpha, rho)
% example
% plot(iwtrans(decfw(wex(3,:)(1:1024),K(1:1024),'d6',4,'wien',[10 400 50 50 50], 4),'d6',4))
% testfwd(3,'d6',4,[10 400 50 50 50],10)
[w, ratiothr] = decfw(wex(anindex,:)(1:1024),K(1:1024),type,p,'wien',alpha, rho);
