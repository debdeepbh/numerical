%% Inputs one signal
%% outputs the quantities related to only the jth wavelet level
function [mean_coeffs_thresholded, ratiounthres, threshold_val]  = multi_fw_one_level(fwax, faximp, fori, noiseax, scaling_j, j)

%tic

global N
global M
global fB
global type
global p
global rho
global method


%% placeholder for all lambdas
lambda = zeros(1, N);


%We shall work for every wavelet level separately
%for j = 1:p+1
    % basis elements
	%%fPhi = fft(B(j,:));
    fPhi = fB(j,:);
    
    % Calulation the avg of all wiener deconvolutions
    % there is not sum and averaging here
        [fdec, lambda(1,:)] = multi_fdecwien(fwax(1,:), faximp(1,:), noiseax(1), fori, scaling_j);
    
    % the inverse Fourier Transform to get signal in time domain
    mean_dec = real(ifft(fdec));
    
    %% Fine till here
    
    % wavelet transfrom
    mean_wdec = wtrans(mean_dec, type, p);
    
    % Choose the coefficients in the jth level
    mean_coeffs = coeff(mean_wdec, p, j);
    
    
    %% Find the threshold
    mean_leaked =  (abs(lambda(1,:)).^2) ./ (abs(faximp(1,:)).^2) .* (noiseax(1).^2);
    % write the sum over frequencies as a dot product
    sigma_sq = dot(mean_leaked, abs(fPhi).^2) / (N);
    threshold_val = sqrt(sigma_sq);
    
    % Apply hard thresholding
    cutoff_vector = (abs(mean_coeffs) > threshold_val);
    mean_coeffs_thresholded = mean_coeffs .* cutoff_vector;
    
    ratiounthres = sum(cutoff_vector)/length(mean_coeffs);
    
    figure;
    stem(mean_coeffs, 'g.'); jupyplot('before')
    figure; 
    stem(mean_coeffs_thresholded,'r.'); jupyplot('after')
    
%printf("unthresholded: %d \t total: %d \t ratio_unthresholded: %f\n", sum(cutoff_vector), length(mean_coeffs), sum(cutoff_vector)/length(mean_coeffs));
    
    
    % place the final coefficients in place
    %[starting, ending] = getindex(N,p,j);
    %w(starting:ending) = mean_coeffs_thres;
%end
        

%toc