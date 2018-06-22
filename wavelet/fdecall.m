% Wiener deconvolution in the frequency domain
function [fw, mult] = fdecall(fsig, fimp, method, scaling)

%sigma = 5;
% Wow, this really works!!
% the noisefile contains: csvwrite('noise_intp.csv',noiseund(320:320+1023))
% where noiseund is from getdata
sigma = csvread('noise_intp.csv');
sigma = abs(fft(sigma));
%sigma = sigma/length(sigma);

%% For testing testy
sigma = abs(fft(randn([1 1024])*5));
% Default option
%sigma = 5;

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
	mult = hsq ./( hsq + scaling*N*(sigma.^2)./(abs(fori).^2));
case 'wien2'
	% definition
	hsq = abs(fimp).^2;

	% model of the original signal
	fori = fsig;
	%fori = fimp;

	% construct the multiplier
	mult = abs(fimp).*hsq ./( hsq + scaling*N*(sigma.^2)./(abs(fori).^2));
otherwise
	disp('wrong method');
end

% multiply by the scaling and take inverse fourier transform
fw = wfft.*mult;


%ww = fw;
%%%%%%%%%%%%%%%%%%%%%%%%
%%% fourier smooth cutoff mutiplier
%mult = zeros(1,length(ww))+1;
%%mult(250:length(ww)-250) = 0;
%mult(length(ww)/8:length(ww)-length(ww)/8) = 0;
%ww = fft(ww).*mult;
%fw = ww;

