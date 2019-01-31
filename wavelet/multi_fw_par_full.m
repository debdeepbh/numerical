%% works well
%% outputs the quantities related to only the jth wavelet level
function w  = multi_fw_par_full(fwax, faximp, fori, noiseax, scaling, j)

%tic

global N
global M
global fB
global type
global p
global rho
global method

w = zeros(1,N);

for j=1:p+1
    [mean_coeffs_thres, ratiounthres, thres_val]  = multi_fw_par_level(fwax, faximp, fori, noiseax, scaling(j), j);
    
    % place the final coefficients in place
    [starting, ending] = getindex(N,p,j);
    w(starting:ending) = mean_coeffs_thres;
    
    ratiounthres
end

%toc