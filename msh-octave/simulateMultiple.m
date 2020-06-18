function [NbdArr_out, u0_multi] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr_multi, Vol_multi, nbd_Vol_multi, extforce_multi, normal_stiffness, contact_radius, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, timesteps);

% break bonds or not
break_bonds = 0;

close all

% get the material properties
%[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta);

[total_nodes, temp, temp] = size(NbdArr_multi);

% pairwise contacts, n choose 2
contact_indices = nchoosek(1:total_particles, 2);
[total_contacts, temp] = size(contact_indices);

% % time integration updaters
% uold_multi = zeros(total_nodes,2, total_particles);
% uolddot_multi = zeros(total_nodes,2, total_particles);
% uolddotdot_multi = zeros(total_nodes,2, total_particles);



imgcounter = 1;
%dt = 25e-9;
%dt = 25e-8;
dt = 0.02/timesteps*2;	% peridem
 
%f = figure('visible','off');
f = figure('visible','on');
tic

for t = 1:timesteps

    u0_multi = uold_multi + dt * uolddot_multi + dt * dt * 0.5 * uolddotdot_multi;

    % placeholder for contact forces on nodes due to neighboring particles
    peridynamic_force_multi = zeros(size(u0_multi));
    contact_force_multi = zeros(size(u0_multi));

%% for debugging, only one particle
    for i=1:1
    %for i=1:total_particles

	% internal force - peridynamic
	[peri_force, stretch] = peridynamic_force(u0_multi(:,:,i), NbdArr_multi(:,:,i), nbd_Vol_multi(:,:,i), xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), cnot) ;

	peridynamic_force_multi(:,:,i) = peri_force;
	stretch_multi(:,:,i) = stretch;



	% contact forces
	if total_particles > 1

%% check: for contact forces, work with current position, instead of deformation u
	    CurrPos_multi = u0_multi + Pos_multi;


	    % information about the central particle
	    CurrPos_center = CurrPos_multi(:,:,i);

	    %% aim to reduce this count, since not all n choose 2 are touching always
	    %for k = 1:total_contacts

	    % placeholder for total contact force on the nodes of i-th body
	    cf_contrib_i = zeros(total_nodes, 2);


	    for j = nonzeros( (1:total_particles).*(1:total_particles ~= i) )

		    % information about the neighboring particle
		    CurrPos_neighbor = CurrPos_multi(:,:,j);
		    Vol_neighbor = Vol_multi(:,:,j);

		    % nodes in the neighboring body that are contact_radius-close to the nodes in the central body
		    [contact_NbdArr] = gen_NbdArr(CurrPos_center, CurrPos_neighbor, contact_radius, 0);
		    % mask for contact neighbors
		    mask = ~~contact_NbdArr;

		    % volume of the contact neighbors
		    cnbd_Vol_neighbor = Vol_neighbor(contact_NbdArr + ~contact_NbdArr) .* mask;

		    CurrPos_neighbor_1 = CurrPos_neighbor(:,1);
		    CurrPos_neighbor_2 = CurrPos_neighbor(:,1);
		    CurrPos_center_1 = CurrPos_center(:,1);
		    CurrPos_center_2 = CurrPos_center(:,1);

		    % direction vectors from the nodes in the central body to the delta-close nodes in the neighboring body
		    direction_1 = (CurrPos_neighbor_1(contact_NbdArr + ~contact_NbdArr) - CurrPos_center_1) .* mask;
		    direction_2 = (CurrPos_neighbor_2(contact_NbdArr + ~contact_NbdArr) - CurrPos_center_2) .* mask;

		    direction_norm = sqrt(direction_1.^2  + direction_2.^2);

		    % in the negative direction, i.e. from the neighbors to the center
		    direction_unit_1 = - direction_1 ./ (direction_norm + ~direction_norm) .* mask;
		    direction_unit_2 = - direction_1 ./ (direction_norm + ~direction_norm) .* mask;


		    % positive contact radius
		    cont_rad_contrib = (contact_radius - direction_norm);
		    positive_cont_rad_contrib = cont_rad_contrib .* (cont_rad_contrib > 0);

		    %% Check this formula to see if the volume is properly accounted for
		    % multiplied by the volumes of contact neighbors
		    contact_force_1 = normal_stiffness .* positive_cont_rad_contrib .* cnbd_Vol_neighbor .* direction_unit_1;
		    contact_force_2 = normal_stiffness .* positive_cont_rad_contrib .* cnbd_Vol_neighbor .* direction_unit_2;


		    total_contact_force = [sum(contact_force_1, 2) sum(contact_force_2,2)];  
		    %% add this force to total internal force
		    cf_contrib_i = cf_contrib_i + total_contact_force;

	    end	% endfor loop in j

	    contact_force_multi(:,:, i) = cf_contrib_i;

	end	%endif

    end %endfor loop in i
    % add all the forces
    total_force_multi = peridynamic_force_multi + contact_force_multi;

    % accelaration
    u0dotdot_multi = (1/rho) .*( total_force_multi + extforce_multi) ;

    % velocity
    u0dot_multi = uolddot_multi + dt * 0.5 * uolddotdot_multi + dt * 0.5 * u0dotdot_multi;

    % update the time loop
    uold_multi = u0_multi;
    uolddot_multi = u0dot_multi;
    uolddotdot_multi = u0dotdot_multi;

    if break_bonds == 1
%% edit to accommodate multiple particles, introduce NbdArr_multi
	    % Bond breaking criteria
	    NbdArr_multi = ((stretch_multi - snot ) < 0).* NbdArr_multi;
    end

    if mod(t, 100) == 0
		t
		    i
		    j
		    disp('total neighbor')
		    sum(sum(mask))

		disp 'contact force contribution'
		max(max(contact_force_multi))
		disp 'peridynamic force contribution'
		max(max(peridynamic_force_multi))


		savenewpos2_multi(total_particles, u0_multi, Pos_multi, 3, imgcounter, f, 'multi_', contact_radius)
		imgcounter = imgcounter + 1;
    end

end

toc

% output
NbdArr_out_multi = NbdArr_multi;




