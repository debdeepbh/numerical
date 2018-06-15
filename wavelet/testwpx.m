% test the values of forward
function testwpx = testwpx(type, p, sigma,  alpha, rho, method)
% example
% plot(iwtrans(decfw(wex(3,:)(1:1024),K(1:1024),'d6',4,'wien',[10 400 50 50 50], 4),'d6',4))
% testfwd(3,'d6',4,[10 400 50 50 50],10)
% optimal match at testfwd(3,'d6',4,[110 100 144 432 320],2)
% Strategy: first get the optimal values using rho = 1
% then adjust rho looking at the wavelet transform
% testspike(testy,'d6',4,5,[0.35 0.35 0.3 0.65 0.7],[2 2 1.5 1 1],'soft')
% For anita data: testspike(wex(3,:), 'd10', 5,sqrt(8.5),[0.3 0.3 0.36 0.64 0.49 4],[4 4 4 1 1 1],'hard')
% bumpless noiseless barebone: testspike(wex(5,:), 'd10', 5,sqrt(8.5),[0.3 0.3 0.36 0.64 0.49 4],[4 4 4 1 3 3],'hard')

close all;
getdata;

for i=1:8
	z = wpx(i,:);
	[w, ratiothr] = wienforwd(z(1:1024),K(1:1024),type,p,sigma,alpha, rho,method);
	figure(1);
	subplot(4,2,i);
	ww = iwtrans(w,type,p);
	plot(ww);
	xlim([1 1024]);
	title(num2str(getsnr(ww)));
	figure(2);
	subplot(4,2,i);
	wf = decall(z,K,'wien',1)(1:1024);
	plot(wf);
	%print -dpng '-S800,1024' wienout.png
	title(num2str(getsnr(wf)));
	xlim([1 1024]);

	%figure(3);
	%subplot(4,2,i);
	%plot(decall(z,K,'allp',1)(1:1024) - iwtrans(w,type,p));
	%title(num2str(getsnr(ww)));
	%xlim([1 1024]);

	%figure(4);
	%subplot(4,2,i);
	%plot(decall(z,K,'naive',1)(1:1024) - iwtrans(w,type,p));
	%title(num2str(getsnr(ww)));
	%xlim([1 1024]);
end

