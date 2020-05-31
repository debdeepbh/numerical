function [NbdArr_out, u0, u_norm, disp] = simulate(Pos, NbdArr, NbdVol, extforce, delta, xi_1, xi_2, xi_norm)

close all

% break bonds or not
break_bonds = 0;

% get the material properties
[rho, Gnot, E, nu, snot, cnot] = material_properties(delta);


total_nodes = length(NbdArr);
% time integration updaters
uold = zeros(total_nodes,2);
uolddot = zeros(total_nodes,2);
uolddotdot = zeros(total_nodes,2);

imgcounter = 1;
dt = 25e-9;

%f = figure('visible','off');
f = figure('visible','on');
tic
for t = 1:1200
    [u0, u0dot, u0dotdot, stretch] = time_update(dt, NbdArr, uold, uolddot, uolddotdot, extforce, NbdVol, rho, cnot, xi_1, xi_2, xi_norm, 'velocity_verlet');

    % update time integration loop
    uold = u0;
    uolddot = u0dot;
    uolddotdot = u0dotdot;
    
    % Bond breaking criteria
    if break_bonds == 1
	    NbdArr = ((stretch - snot ) < 0).* NbdArr;
    end

    % plotting
    if mod(t, 50) == 0
	t
	savenewpos2(u0, Pos, 3, imgcounter, f, 'test')
    	imgcounter = imgcounter + 1;
    end


end

toc

% output
NbdArr_out = NbdArr;

%drawmesh(out, dx, dy, nx, ny);



