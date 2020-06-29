function [delta, rho, Gnot, E, nu, snot, cnot, bulk_modulus] = material_properties(delta, material) 
% Function description
% Input:
%	: 
%
% Output:
%	sigma:
%	rho:
%	Gnot:
%	E:
%	nu:
%	snot:
%	cnot:


%initial condition
%sigma=12e6/dx;


switch material
    case 'sodalime'
	% Original
	E = 72e9;
	rho=2440;
	nu = 0.22;       % commnet on paper: the effective poisson ration is 1/3

	cnot = 6*E/( pi * (delta^3) * (1 - nu));

	Gnot = 135;
	snot = sqrt(4 * pi * Gnot /(9*E*delta));

	bulk_modulus = E/ (3 * ( 1 - 2 *nu)); 

    case 'peridem'
	%% peridem simulation
  % K 
  bulk_modulus = 2.e+09;
  % G
  shear_modulus = 1.33e+09;
  rho=1200	% peridem

nu = (3 * bulk_modulus - 2 * shear_modulus) / ( 2 * ( 3 * bulk_modulus + shear_modulus));
% nu = 0.2278
E = 9 * bulk_modulus * shear_modulus / ( 9 * bulk_modulus + shear_modulus);
% E = 1.23e9

	%lambda = E * nu / ((1 + nu)*( 1 - 2 * nu));
	%cnot = 24 * lambda/ (pi * delta^3);


	cnot = 24 * E /( (1 - nu) * pi * delta^3);

	Gnot = 135;
	snot = sqrt(4 * pi * Gnot /(9*E*delta));
		
    otherwise
	
end




