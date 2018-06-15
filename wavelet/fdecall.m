% Wiener deconvolution in the frequency domain
function [fw, mult] = fdecall(fsig, fimp, method, scaling)

sigma = 5;

N = length(fsig);
if ~exist('scaling')
	scaling = 1;
end

% assume that the fsig and fimp are of the same length

% naive deconvolution in the fourier domain
wfft = fsig./fimp;

switch method
case 'naive'
	mult = 1;
case 'allp'
	mult = abs(fimp);
case 'wien'
	% definition
	hsq = abs(fimp).^2;

	% model of the original signal
	fori = fsig;
	%fori = fimp;

	% construct the multiplier
	mult = hsq ./( hsq + scaling*N*(sigma^2)./(abs(fori).^2));
case 'wien2'
	% definition
	hsq = abs(fimp).^2;

	% model of the original signal
	fori = fsig;
	%fori = fimp;

	% construct the multiplier
	mult = abs(fimp).*hsq ./( hsq + scaling*N*(sigma^2)./(abs(fori).^2));
otherwise
	disp('wrong method');
end

% multiply by the scaling and take inverse fourier transform
fw = wfft.*mult;

