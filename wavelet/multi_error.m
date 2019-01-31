%% Failed to design an experiment with so many dependent and independent variables
%% Moving to level-independent
%% compute the error for tikhonov regularization


% set the parameters
type = 'meyer';
p = 3;
rule = 'hard';
rho  = 1;

M = 15;

multi_data_same

N = length(wax);
fori = fft(testyori);
w = wtrans(testyori, type, p);

wFinal = zeros(1,N);
wFinalIndi = zeros(M,N);

% iterations
tauList = 0.001:0.01:0.2;
%tauList = 0.5;
totalIter = length(tauList);

MSE = zeros(p,totalIter);
actualErrLevel = zeros(1,totalIter);

for tauIndex = 1:totalIter 
%for tauIndex = 1 

for j=1:p+1
%for j = 4	% wavelet level



B = getbasismat(type, p, N);



fPhi = fft(B(j,:));
wcoeffs = coeff(w, p, j);
Nj = length(wcoeffs);

	tau = tauList(tauIndex);
	freg = zeros(M,N);
	sigmalHolder = 0;
	decHolder = zeros(M,N);
	for i=1:M

		sig = wax(i,:);
		fsig = fft(sig);
		imp = aximp(i,:);
		fimp = fft(imp);
		sigma = std(noiseax(i,:));
		%fori = length(wax)* (sigma^2)/alpha; 

		% perform Tikhonov decovolution with alpha_j*N*sigma^2/|X|^2 = tau;
		[fdec, mult] = multi_fdecwien(fsig, fimp, sigma, fori, tau);
		decHolder(i,:) = fdec;
		%[fdec, mult] = multi_fdecwien(fsig, fimp, sigma, fori, scaling(j));

		% save the fourier regularization parameter
		freg(i,:) = mult;

		sigmalHolder = sigmalHolder + sigma^2 * dot((abs(mult).^2) ./ (abs(fimp).^2), abs(fPhi).^2);

		% sigmal for individual ForWarD purposes
		sigmalIndi(i) = sigma^2/N*dot(abs(mult).^2, (abs(fPhi).^2) ./ (abs(fimp).^2));

	end

	%% using level-dependent expression
	MSE1 = Nj * dot((abs(fori).* abs(fPhi)).^2, (1 - mean(freg)).^2)/N;
	%% change it here to compute the mean error


	%%% compute sigmal, level dependent effective leaked noise
	sigmal = sqrt(sigmalHolder /N/(M^2));

	MSE2 = sum(min(abs(wcoeffs).^2, sigmal^2));

	%z_par = multi_fpar(wax, aximp, testyori, type, p,  noiseax, alpha_par, rho, rule);
	MSE(j,tauIndex) = MSE1 + MSE2;


	%%%% actual ForWaRD
	%% averaged %%
	favgDec = mean(decHolder); 
	wLevel = coeff(wtrans(real(ifft(favgDec)), type, p), p, j); 


	% applythresshold and zero out irrerelevant coeffs
	%% hard thresholding
	wLevelThres = (abs(wLevel) >= sigmal).*wLevel;

	% storing the thresholded vector in appropriate location of output
	[indF, indL] = getindex(N, p, j);
	wFinal(indF:indL) = wLevelThres;

	%%%%%% individual Forward%%
	wLevelIndi = zeros(1, Nj);
	wLevelThresIndi = zeros(M, Nj);
	for i = 1:M
		% getting the j-th level coefficients for wtrans
		wLevelIndi(1,:) = coeff(wtrans(real(ifft(favgDec)), type, p), p, j); 
		%% hard thresholding
		wLevelThresIndi(i,:)= (abs(wLevelIndi) >= sigmalIndi(i)).*wLevelIndi;
		
		% storing the thresholded vector in appropriate location of output
		wFinalIndi(i,indF:indL) = wLevelThresIndi(i,:);
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	end


	
end



	%% relative error of the wavelet coefficients
	% j-th level error for geenralized Forward
	actualErrLevel(tauIndex) = rele(wLevelThres, wcoeffs);
	% j-th level error of individual 
	for i=1:M
		actualErrLevelIndi(i,tauIndex) = rele(wLevelThresIndi(i,:), wcoeffs);
	end
	% j-th level error for the mean of individually deconvolved
	actualErrLevelIndiMean(tauIndex) = rele(mean(wLevelThresIndi), wcoeffs);




end

for j=1:p+1
	subplot(p+1, 1, j)
	plot(tauList, MSE(j,:))
end
%plot(tauList, actualErrLevel, tauList, actualErrLevelIndiMean);
%legend('GForWd', 'IndiMean')


%%% plot error %%%%%%%%%%%
%plot(tauList, log(MSE));
%hold on

%plot(tauList, MSE);
%hold off

%%%%%%%% plot output pictures %%%%
%% plot generalized forward
%figure 1
%plot(iwtrans(wFinal, type, p))
%hold on
%
%% plot individually deconvolved ones
%figure 2;
%hold on
%for i =1:M
%	plot(iwtrans(wFinalIndi(i,:), type, p));
%end
%hold off
%
%% plot avg of individually deconvolved
%figure 1
%plot(iwtrans(mean(wFinalIndi), type, p))
%hold off
%


