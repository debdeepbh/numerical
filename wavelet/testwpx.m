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

% modifying K here
%K = real(ifft(fft(K)./abs(fft(K))));
%K = compress(K, 'd10', 4, 50);

for i=1:8
	z = wpx(i,:);
	z = z(1:1024);
	[w, ratiothr, thrvec] = wienforwd(z,K,type,p,sigma,alpha, rho,method);
	
	switch type
	case 'bior13d'
		type = 'bior13r'
	case 'bior33d'
		type = 'bior33r'
	end
	ww = iwtrans(w, type,p);



% take projection 
% same as reconstruction upto level 2, which gives 1250MHz
%%%%%%%%%%%%%%%%%%%%%%%
ww = proj(ww, type,2);
%%% fourier cutoff mutiplier (windowing)
%mult = zeros(1,length(ww))+1;
%%mult(250:length(ww)-250) = 0;
%mult(length(ww)/8:length(ww)-length(ww)/8) = 0;
%ww = real(ifft(fft(ww).*mult));


	figure(1);
	subplot(4,2,i);
	plotsnr(ww);
	xlim([1 1024]);

	figure(2);
	subplot(4,2,i);
	wf = decwien(z,K,sigma,1);
	plotsnr(wf);
	xlim([1 1024]);

	figure(3);
	subplot(4,2,i);
	plotpanita(ww);
	

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

