% Wiener deconvolution in the frequency domain
function [fw, mult] = fwiendec(fsig, fimp, sigma, scaling)

N = length(fsig);
if ~exist('scaling')
	scaling = 1;
end

% assume that the fsig and fimp are of the same length

% naive deconvolution in the fourier domain
wfft = fsig./fimp;
% definition
hsq = abs(fimp).^2;

% model of the original signal
fori = fsig;
%fori = fimp;

% construct the multiplier
mult = hsq ./( hsq + scaling*N*(sigma.^2)./(abs(fori).^2));
% multiply by the scaling and take inverse fourier transform
fw = wfft.*mult;

