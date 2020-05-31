function [u0, u0dot, u0dotdot, stretch] = update_timeint(uold, uolddot, uolddotdot, dt, NbdArr, Vol, xi_1, xi_2, xi_norm, extforce,  cnot, rho) 
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

% loop starts
u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;

disp_1 = u0(:,1);
disp_2 = u0(:,2);

% Displacements of the centers
%u_1 = disp_1(longvec);
%u_2 = disp_2(longvec);

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

%unitvec_1 =  eta_p_xi_1 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* (~~eta_p_xi_norm);
%unitvec_2 =  eta_p_xi_2 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* (~~eta_p_xi_norm);

Force1 = cnot.* stretch .* unitvec_1;
Force2 = cnot.* stretch .* unitvec_2;

%% mask it
Force1 = Force1 .* restrict;
Force2 = Force2 .* restrict;

totalintforce = [sum(Force1.*restrict, 2) sum(Force2.*restrict,2)];    %sum(A) ad totalintforce2 = sum(Force2, 2);

u0dotdot = (1/rho) .*( totalintforce .* Vol + extforce) ;

%calculate u0dot
u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;
