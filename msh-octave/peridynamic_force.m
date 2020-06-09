function [totalintforce, stretch] = peridynamic_force(u0, NbdArr, nbd_Vol, xi_1, xi_2, xi_norm, cnot) 
% Function description
% Input:
%	uold: 
%	uolddot: 
%	uolddotdot: 
%	: 
%
% Output:
%	out:


restrict = (NbdArr > 0);

disp_1 = u0(:,1);
disp_2 = u0(:,2);

%% Displacements of the neighbors
u_1_p = disp_1(NbdArr + (~NbdArr));
u_2_p = disp_2(NbdArr + (~NbdArr));

eta_1 = u_1_p - disp_1;
eta_2 = u_2_p - disp_2;

eta_p_xi_1 = eta_1 + xi_1;
eta_p_xi_2 = eta_2 + xi_2;

eta_p_xi_norm = sqrt(eta_p_xi_1.^2 + eta_p_xi_2.^2);

% to avoid dividing by zero 
stretch = (eta_p_xi_norm - xi_norm)./(xi_norm + (~xi_norm));

unitvec_1 =  eta_p_xi_1 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* restrict;
unitvec_2 =  eta_p_xi_2 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* restrict;

Force1 = cnot.* stretch .* unitvec_1;
Force2 = cnot.* stretch .* unitvec_2;

% nbd_Vol is already masked
totalintforce = [sum(Force1.* nbd_Vol, 2) sum(Force2.* nbd_Vol,2)];  

