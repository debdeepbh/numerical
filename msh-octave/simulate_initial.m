function [NbdArr_out, u0] = simulate_initial(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, extforce, delta, xi_1, xi_2, xi_norm, timesteps)

% break bonds or not
break_bonds = 0;

close all

% get the material properties
[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta);

total_nodes = length(NbdArr);
% % time integration updaters
%uold = zeros(total_nodes,2);
%uolddot = zeros(total_nodes,2);
%uolddotdot = zeros(total_nodes,2);

imgcounter = 1;
%dt = 25e-9;
dt = 25e-8;

f = figure('visible','off');
%f = figure('visible','on');
tic

for t = 1:timesteps
    %for t = 1:2800

    [u0, u0dot, u0dotdot, stretch] = update_timeint(uold, uolddot, uolddotdot, dt, NbdArr, Vol, xi_1, xi_2, xi_norm, extforce, cnot, rho);

    uold = u0;
    uolddot = u0dot;
    uolddotdot = u0dotdot;

    if break_bonds == 1
	    % Bond breaking criteria
	    NbdArr = ((stretch - snot ) < 0).* NbdArr;
    end

    if mod(t, 50) == 0
		t
		savenewpos2(u0, Pos, 3, imgcounter, f, 'test')
		imgcounter = imgcounter + 1;
    end

end

toc

% output
NbdArr_out = NbdArr;




