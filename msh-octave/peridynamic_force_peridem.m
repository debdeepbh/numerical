function [totalintforce, stretch ] = peridynamic_force_peridem(CurrPos, NbdArr, nbd_Vol, xi_1, xi_2, xi_norm, cnot, delta, i, t) 
% Function description
% Input:
%	uold: 
%	uolddot: 
%	uolddotdot: 
%	: 
%
% Output:
%	out:

% cnot = 24 * E/ (1 - nu)
% J = 1/ ( 1 - xi_norm/delta)


restrict = (NbdArr > 0);

%%%%%%%%%%5
CurrPos_1 = CurrPos(:,1);
CurrPos_2 = CurrPos(:,2);

eta_p_xi_1 = CurrPos_1(NbdArr + ~NbdArr)  - CurrPos_1;
eta_p_xi_2 = CurrPos_2(NbdArr + ~NbdArr)  - CurrPos_2;

eta_p_xi_1 = eta_p_xi_1 .* restrict;
eta_p_xi_2 = eta_p_xi_2 .* restrict;

eta_p_xi_norm = sqrt(eta_p_xi_1.^2 + eta_p_xi_2.^2);

    if ((i==2) && (t == 13))
	i
	t
	max(max(CurrPos_2))
    end

% to avoid dividing by zero 
stretch = (eta_p_xi_norm - xi_norm)./(xi_norm + (~xi_norm)) .* (~~xi_norm);
% zero out negative stretch
%% debug
stretch = stretch .* (stretch > 0);

unitvec_1 =  eta_p_xi_1 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* (~~eta_p_xi_norm);
unitvec_2 =  eta_p_xi_2 ./ (eta_p_xi_norm + ~eta_p_xi_norm) .* (~~eta_p_xi_norm);

% influence function J = 1/(1-r)
J_denom = 1 - xi_norm/delta;
J = 1 ./ (J_denom + ~J_denom) .* (~~J_denom);


Force1 = cnot .* J .* stretch .* unitvec_1;
Force2 = cnot .* J .* stretch .* unitvec_2;


% nbd_Vol is already masked
totalintforce = [sum(Force1.* nbd_Vol, 2) sum(Force2.* nbd_Vol,2)];  

