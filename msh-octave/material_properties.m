function [delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta) 
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

%experiment = 'original'
experiment = 'peridem'

switch experiment
    case 'original'
	% Original
	E = 72e9;
	rho=2440;
	nu = 0.25;       % commnet on paper: the effective poisson ration is 1/3

	cnot = 6*E/( pi * (delta^3) * (1 - nu));

	Gnot = 135;
	snot = sqrt(4 * pi * Gnot /(9*E*delta));

    case 'peridem'
	%% peridem simulation
	E = 3.24e9;	% peridem
	rho=1200	% peridem
	nu = 0.22;	% peridem

	lambda = E * nu / ((1 + nu)*( 1 - 2 * nu));
	cnot = 24 * lambda/ (pi * delta^3);


	Gnot = 135;
	snot = sqrt(4 * pi * Gnot /(9*E*delta));
		
    otherwise
	
end




