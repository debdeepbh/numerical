% different deconvolution methods
function w = decall(signal, impulse, method, scaling)
% scaling is optional
% var is also tau for tikhonov

% figure out optimal padding
% pad impulse by zero to make them of the same length
% assuming here that impulse is of smaller length than signal
N = length(signal);
paddedimpulse = zeros(1,N);
%starting = floor((length(signal)-length(impulse))/2);
starting=1;
paddedimpulse(starting:starting+length(impulse)-1) = impulse;

    
%     sigpower = abs(fft(signal)).^2;
%     % fake noise power value from the signal
%     noisepower = (norm(signal(801:1000),2)^2)/sqrt(200)*length(signal);
%     snrval = (sigpower)./(noisepower);
    
    % constant snrval
    % rough extimate (possibly wrong) from pspectrum
    % But above 0.5 it does not even matter
 %snrval = 10;% abs(fft(signal)).^2/2; 
 %snrval = 2; %csvread('snrfromobs.csv') ;   


% naive deconvolution in the fourier domain
fimp = fft(paddedimpulse);
fsig = fft(signal);
wfft = fsig./fimp;

% construct the fourier multiplier based on method
switch method
case 'none'
	% do nothing
	mult = 1;

case 'allp'	% allpass
	mult = abs(fimp);

case 'tikh'	% tikhnov, constant tau
    hsq = abs(fimp).^2;
    % in fact, this var is actually tau
    % but this can be set using var, abuse of notation
    mult = hsq ./( hsq + var);

case 'wien'	% wiener, following Neelamani 2004
   
    var = 8.5;	% get this from estimating the noise

    hsq = abs(fimp).^2;
    mult = hsq ./( hsq + scaling*N*(var^2)./(fsig.^2));
    
otherwise
	printf('Wrong method');
end


% multiply by the scaling and take inverse fourier transform
w = real(ifft(wfft.*mult));

