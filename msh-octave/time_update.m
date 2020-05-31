function [u0, u0dot, u0dotdot, stretch] = time_update(dt, NbdArr, uold, uolddot, uolddotdot, extforce, NbdVol, rho, cnot, xi_1, xi_2, xi_norm, time_integration_scheme) 
% Performs one step of time integration. Given the displacement, velocity, and acceleration at the current time step, returns the same vectors for the next time step.
% n = # nodes, l = #max neighbor
% Input:
%	uold: old displacement vector, nx2
%	uolddot: old velocity vector, nx2 
%	uolddotdot: old acceleration vector, nx2
%	extforce: external force being applied to the nodes, nx2
%	NbdVol: nodal volume of the neighbors, nxl 
%	rho: density of the material (or of the nodes), scalar (or nx1)
%	xi_1, xi_2, xi_norm: vectors associated with xi = (x' - x), nxl each
%	time_integration_scheme: scheme to perform the time integration
%
% Output:
%	u0: new displacement vector, nx2
%	u0dot: new velocity vector, nx2 
%	u0dotdot: new acceleration vector, nx2
%	stretch: stretch of the neighbors, nxl


%restrict = (Nbd > 0);
restrict = ~~NbdArr;

switch time_integration_scheme
case 'velocity_verlet'
	% loop starts
	u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;


	% dispacement vector
	disp_1 = u0(:,1);
	disp_2 = u0(:,2);

	%% Displacements of the neighbors
	u_1_p = disp_1(NbdArr + (~NbdArr));
	u_2_p = disp_2(NbdArr + (~NbdArr));

	% displacements of the neighbors relative to the center
	eta_1 = u_1_p - disp_1;
	eta_2 = u_2_p - disp_2;


	eta_plus_xi_1 = eta_1 + xi_1;
	eta_plus_xi_2 = eta_2 + xi_2;

	eta_plus_xi_norm = sqrt(eta_plus_xi_1.^2 + eta_plus_xi_2.^2);

	% to avoid dividing by zero 
	stretch = (eta_plus_xi_norm - xi_norm)./(xi_norm + (~xi_norm));

	unitvec_1 =  eta_plus_xi_1 ./ (eta_plus_xi_norm + ~eta_plus_xi_norm);
	unitvec_2 =  eta_plus_xi_2 ./ (eta_plus_xi_norm + ~eta_plus_xi_norm);

	Force1 = cnot.* stretch .* unitvec_1;
	Force2 = cnot.* stretch .* unitvec_2;

	%% mask it
	Force1 = Force1 .* restrict;
	Force2 = Force2 .* restrict;


	% total internal force is obtained by adding the forces associated with the neighbors
	% Summing over forces along with the nodal volumes, automatically masked
	totalintforce = [sum(Force1 .* Vol, 2) sum(Force2 .* Vol, 2)];   

	u0dotdot = (1 ./rho) .* ( totalintforce + extforce) ;

	u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;

otherwise
	disp 'Wrong scheme specified'
end
