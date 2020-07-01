function [NbdArr_out, u0] = test_int(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, nbd_Vol, extforce, delta, xi_1, xi_2, xi_norm, timesteps, break_bonds, material, dt)
% Function description
% Input:
%	uold: 
%	udotold: 
%	udotdotold: 
%
% Output:
%	out:
imgcounter = 1;

scheme = 'verlet'

%% for gravity test
%dt = 1e-3;

%f = figure('visible','off');
f = figure('visible','on');

[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta, material);

totalintforce = zeros(size(extforce));

for t = 1:timesteps
    %for t = 1:2800
    t

    switch scheme
        case 'verlet'
	    u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;

	    u0dotdot = (1/rho) .*( totalintforce  + extforce) ;

	    %max(max(u0dotdot))

	    % velocity
	    u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;

	    %update loop
	    uold = u0;
	    uolddot = u0dot;
	    uolddotdot = u0dotdot;


	case 'eular'
    	
        otherwise
    	
    end

    if mod(t, 50) == 0
		t
		time = dt * t * 1e3 ; % milli second
	    % % Check for constant acceleration
	    %yval = u0(:,2);
	    %fprintf('Is: %f\n', (min(yval) + max(yval))/2);
	    %fprintf('Should be: %f\n', 1/2 * (extforce(1,2)/rho) * (t*dt)^2);

		scaling = 1;
		savenewpos2(u0, Pos, 3, imgcounter, f, 'single_', time, scaling)
		imgcounter = imgcounter + 1;
    end

end

NbdArr_out = NbdArr;
