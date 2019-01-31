%% works well
%% outputs the quantities related to only the jth wavelet level
function [mean_coeffs_thres, ratiounthres, threshold_val]  = multi_fw_par_level(fwax, faximp, fori, noiseax, scaling_j, j)

%tic

global N
global M
global fB
global type
global p
global rho
global method


%% placeholder for all lambdas
lambda = zeros(M, N);


%We shall work for every wavelet level separately
%for j = 1:p+1
    % basis elements
	%%fPhi = fft(B(j,:));
    fPhi = fB(j,:);
    
    % Calulation the avg of all wiener deconvolutions
    mean_fded = sumi_fdec = zeros(1,N);
    for i=1:M
        % Why are we saving lambda for all i and j?
        [fdec, lambda(i,:)] = multi_fdecwien(fwax(i,:), faximp(i,:), noiseax(i), fori, scaling_j);
        sumi_fdec = sumi_fdec + fdec;
    end
    mean_fdec = sumi_fdec/M;
    
    % the inverse Fourier Transform to get signal in time domain
    mean_dec = real(ifft(fdec));
    
    %% Fine till here
    
    % wavelet transfrom
    mean_wdec = wtrans(mean_dec, type, p);
    
    % Choose the coefficients in the jth level
    mean_coeffs = coeff(mean_wdec, p, j);
    
    
    %% Find the threshold
    sum_leaked = mean_leaked = zeros(1,N);
    for i=1:M
        sum_leaked = sum_leaked .+ (abs(lambda(i,:)).^2) ./ (abs(faximp(i,:)).^2) .* (noiseax(i).^2);
    end
    mean_leaked = sum_leaked/M;
        
    % write the sum over frequencies as a dot product
    sigma_sq = dot(mean_leaked, abs(fPhi).^2) / (N*M);
    
    threshold_val = sqrt(sigma_sq);
    
    % Apply hard thresholding
    cutoff_vector = (abs(mean_coeffs) > threshold_val);
    mean_coeffs_thres = mean_coeffs .* cutoff_vector;
    
    ratiounthres = sum(cutoff_vector)/length(mean_coeffs);
    
%printf("unthresholded: %d \t total: %d \t ratio_unthresholded: %f\n", sum(cutoff_vector), length(mean_coeffs), sum(cutoff_vector)/length(mean_coeffs));
    
    
    % place the final coefficients in place
    %[starting, ending] = getindex(N,p,j);
    %w(starting:ending) = mean_coeffs_thres;
%end
        

%toc
