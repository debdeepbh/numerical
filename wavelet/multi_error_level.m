function errvals = multi_error_level(deconv, f_wax, f_aximp, f_ori, w_ori, noise_sd, alphavals, level_now)

global N;
global p;

errvals = zeros(1,length(alphavals));

counter = 1;
for alpha=alphavals
    if deconv == 'par'
        [w_level, first_unthresratio, first_thresval] = multi_fw_par_level(f_wax, f_aximp, f_ori, noise_sd, alpha, level_now);
    elseif deconv == 'one'
        [w_level, first_unthresratio, first_thresval] = multi_fw_one_level(f_wax, f_aximp, f_ori, noise_sd, alpha, level_now);
    end
    
    
    % getting the corresponding wavelet coefficients of the original signal
    [starting, ending] = getindex(N,p,level_now);
    w_ori_level = w_ori(starting:ending);
    
    % All are scaled, but variable for signal or level
    %%computing the relative error
    %errvals(counter) = rele(w_level,w_ori_level);
    %% computing the mean error
    %errvals(counter) = norm(w_level-w_ori_level,2)/length(w_level);
    % computing the l^2 error
    errvals(counter) = norm(w_level - w_ori_level,2);
    counter ++;
end

% plots too
plot(alphavals, errvals)