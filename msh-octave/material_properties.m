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


rho=2440;
Gnot = 135;
E = 72e9;
nu = 0.22;       % commnet on paper: the effective poisson ration is 1/3

snot = sqrt(4 * pi * Gnot /(9*E*delta));
cnot = 6*E/( pi * (delta^3) * (1 - nu));
