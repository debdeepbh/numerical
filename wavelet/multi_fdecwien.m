% Wiener deconvolution in the frequency domain
% assume that the fsig and fimp are of the same length
function [fw, mult] = multi_fdecwien(fsig, fimp, noise, fori, scaling)


N = length(fsig);


	% definition
	hsq = abs(fimp).^2;

	% construct the multiplier
	if (length(noise)==1)	% noise is actually the  sd
		Nsigmasq = N*noise^2;
	else			% it is the raw noise
		Nsigmasq = abs(fft(noise)).^2;
	end
	fw = fsig .* conj(fimp) ./( hsq + scaling*(Nsigmasq)./(abs(fori).^2) );
	% we need to output multi as well, for wienforwd.m
	mult = hsq ./( hsq + scaling*(Nsigmasq)./(abs(fori).^2) );
% multiply by the scaling and take inverse fourier transform
%fw = wfft.*mult;
