%% works well
%% Barebone function, input is the fft, output is wavelet transform
function w = multi_fw_par_new(fwax, faximp, fori, noiseax, scaling)

%tic

global N
global M
global fB
global type
global p
global rho
global method

% output
w = zeros(1,N);

%% placeholder for all lambdas
lambda = zeros(M, p+1, N);


%We shall work for every wavelet level separately
for j = 1:p+1
    % basis elements
	%%fPhi = fft(B(j,:));
    fPhi = fB(j,:);
    
    % Calulation the avg of all wiener deconvolutions
    mean_fded = sumi_fdec = zeros(1,N);
    for i=1:M
        % Why are we saving lambda for all i and j?
        [fdec, lambda(i,j,:)] = multi_fdecwien(fwax(i,:), faximp(i,:), noiseax(i), fori, scaling(j));
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
        %% this extra line is for dimensional mismatch
        % apprantly size(lambda(i,j,:)) = [1 1 1024] and not [1 1024]
        % alternatively, squeeze() drops any dimension with value 1
        current_lambda(1,:) = lambda(i,j,:);
        sum_leaked = sum_leaked .+ (abs(current_lambda).^2) ./ (abs(faximp(i,:)).^2) .* (noiseax(i).^2);
    end
    mean_leaked = sum_leaked/M;
        
    % write the sum over frequencies as a dot product
    sigma_sq = dot(mean_leaked, abs(fPhi).^2) / (N*M);
    
    thres = sqrt(sigma_sq);
    
    % Apply hard thresholding
    cutoff_vector = abs(mean_coeffs) > thres; 
    mean_coeffs_thres = mean_coeffs .* cutoff_vector;
    
printf("unthresholded: %d \t total: %d \t ratio_unthresholded: %f\n", sum(cutoff_vector), length(mean_coeffs), sum(cutoff_vector)/length(mean_coeffs));
    
    % place the final coefficients in place
    [starting, ending] = getindex(N,p,j);
    w(starting:ending) = mean_coeffs_thres;
end
        

%toc