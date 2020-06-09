function [NbdArr_out, u0_multi] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr, Vol_multi, extforce_multi, delta, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, timesteps);

% break bonds or not
break_bonds = 0;

close all

% get the material properties
%[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta);

total_nodes = length(NbdArr);


% % time integration updaters
% uold_multi = zeros(total_nodes,2, total_particles);
% uolddot_multi = zeros(total_nodes,2, total_particles);
% uolddotdot_multi = zeros(total_nodes,2, total_particles);

imgcounter = 1;
%dt = 25e-9;
dt = 25e-8;

f = figure('visible','off');
%f = figure('visible','on');
tic

for t = 1:timesteps
    %for t = 1:2800

    for i=1:total_particles
	[u0_multi(:,:,i), u0dot_multi(:,:,i), u0dotdot_multi(:,:,i), stretch_multi(:,:,i)] = update_timeint(uold_multi(:,:,i), uolddot_multi(:,:,i), uolddotdot_multi(:,:,i), dt, NbdArr, Vol_multi(:,:,i), xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), extforce_multi(:,:,i), cnot, rho);
    end

    uold_multi = u0_multi;
    uolddot_multi = u0dot_multi;
    uolddotdot_multi = u0dotdot_multi;

%% edit to multi particle
    if break_bonds == 1
	    % Bond breaking criteria
	    NbdArr = ((stretch - snot ) < 0).* NbdArr;
    end

    if mod(t, 50) == 0
		t
		savenewpos2_multi(total_particles, u0_multi, Pos_multi, 3, imgcounter, f, 'test')
		imgcounter = imgcounter + 1;
    end

end

toc

% output
NbdArr_out = NbdArr;




