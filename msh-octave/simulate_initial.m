function [NbdArr_out, u0] = simulate_initial(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, nbd_Vol, extforce, delta, xi_1, xi_2, xi_norm, timesteps, break_bonds)


close all

% get the material properties
[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta);

total_nodes = length(NbdArr);
% % time integration updaters
%uold = zeros(total_nodes,2);
%uolddot = zeros(total_nodes,2);
%uolddotdot = zeros(total_nodes,2);


imgcounter = 1;
dt = 25e-9;
%dt = 25e-8;

%f = figure('visible','off');
f = figure('visible','on');
tic

for t = 1:timesteps
    %for t = 1:2800

    % loop starts
    u0 = uold + dt * uolddot + dt * dt * 0.5 * uolddotdot;

    [totalintforce, stretch] = peridynamic_force(u0, NbdArr, nbd_Vol, xi_1, xi_2, xi_norm, cnot);

    % accelaration, combine all the forces
% debuggin: include Volume
    %u0dotdot = (1/rho) .*( totalintforce .* Vol + extforce) ;
    u0dotdot = (1/rho) .*( totalintforce  + extforce) ;

    % velocity
    u0dot = uolddot + dt * 0.5 * uolddotdot + dt * 0.5 * u0dotdot;


    uold = u0;
    uolddot = u0dot;
    uolddotdot = u0dotdot;

    if break_bonds == 1
	    % Bond breaking criteria
	    NbdArr = ((stretch - snot ) < 0).* NbdArr;
    end

    if mod(t, 50) == 0
		t
		savenewpos2(u0, Pos, 3, imgcounter, f, 'single_')
		imgcounter = imgcounter + 1;
    end

end

toc

% output
NbdArr_out = NbdArr;




