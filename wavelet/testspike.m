% test the values of forward, ww is the final result, w is wtrans derived from wienforwd
function [ww, w] =  testspike(z, K, type, p, sigma,  alpha, rho, method, post)

% example
% plot(iwtrans(decfw(wex(3,:)(1:1024),K(1:1024),'d6',4,'wien',[10 400 50 50 50], 4),'d6',4))
% testfwd(3,'d6',4,[10 400 50 50 50],10)
% optimal match at testfwd(3,'d6',4,[110 100 144 432 320],2)
% Strategy: first get the optimal values using rho = 1
% then adjust rho looking at the wavelet transform
% testspike(testy,'d6',4,5,[0.35 0.35 0.3 0.65 0.7],[2 2 1.5 1 1],'soft')
% For anita data: testspike(wex(3,:), 'd10', 5,sqrt(8.5),[0.3 0.3 0.36 0.64 0.49 4],[4 4 4 1 1 1],'hard')
% bumpless noiseless barebone: testspike(wex(5,:), 'd10', 5,sqrt(8.5),[0.3 0.3 0.36 0.64 0.49 4],[4 4 4 1 3 3],'hard')

%getdata;
close all;

% modifying  K  here
%K = real(ifft(fft(K)./abs(fft(K))));
%K = compress(K, 'd10', 5, 50); 
% additional smoothing of the data


[w, ratiount, thrvec] = wienforwd(z(1:1024),K(1:1024),type,p,sigma,alpha, rho,method);

plotthr(w,p,thrvec);
%plotcoeffs(w,p);

%subplot(211)
%plot(z)
%title('observed signal')
%xlim([1 length(w)])

switch type
case 'bior13d'
	type = 'bior13r'
case 'bior33d'
	type = 'bior33r'
end
%subplot(212)
ww = iwtrans(w,type,p);

switch post	% post processing
case 'yes'
%% take projection 
%%%%%%%%%%%%%%%%%%%%%%%%
ww = proj(ww, type, 2);	 % for Meyer's, reconstruction at level 2
			% removes frequencies outside 1500MHz
			% But 800-1500 get removed a bit
%%%% fourier smooth cutoff mutiplier
%mult = zeros(1,length(ww))+1;
%mult(250:length(ww)-250) = 0;
%ww = real(ifft(fft(ww).*mult));
end

%% Plot the deconvolved signal
figure;
plot(ww);
titlestr = strcat('deconvolved signal; method:', method, ', wavelet:', type);
title(titlestr)
xlim([1 length(w)])


%% Plot the power of the deconvolved signal 
%figure;
%%plot(log10(abs(fft(ww)).^2));
%%plot((1:length(ww))*5000/512, abs(fft(ww)));
%plotfanita(ww);

% plot after post processing
%plotcoeffs(wtrans(ww, type,p),p);


% print the threshold proportion
ratiount
thrvec
