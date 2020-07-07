function [NbdArr_out_multi, u0_multi, store_location, store_vel_min, store_vel_max  ] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr_multi, Vol_multi, nbd_Vol_multi, extforce_multi, normal_stiffness, contact_radius, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, dt, timesteps, delta, modulo, break_bonds);


% break bonds or not

%save_plot = 0;
save_plot = 1;

%use_influence_function = 0
use_influence_function = 1

close all

[total_nodes, temp, temp] = size(NbdArr_multi);

% pairwise contacts, n choose 2
contact_indices = nchoosek(1:total_particles, 2);
[total_contacts, temp] = size(contact_indices);

imgcounter = 1;
%dt = 25e-9;
%dt = 25e-8;

%f = figure('visible','off');
f = figure('visible','on');
tic

for t = 1:timesteps

    u0_multi = uold_multi + dt * uolddot_multi + dt * dt * 0.5 * uolddotdot_multi;

    % computing current position
    CurrPos_multi = u0_multi + Pos_multi;

    % placeholder for contact forces on nodes due to neighboring particles
    peridynamic_force_multi = zeros(size(u0_multi));
    contact_force_multi = zeros(size(u0_multi));

%% for debugging, only one particle
    %for i=2
    for i=1:total_particles

	% internal force - peridynamic
	[peridynamic_force_multi(:,:,i), stretch_multi(:,:,i)] = peridynamic_force_bypos(CurrPos_multi(:,:,i), NbdArr_multi(:,:,i), nbd_Vol_multi(:,:,i), xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), cnot, delta, use_influence_function) ;

	% contact forces
	if total_particles > 1

	    % information about the central particle
	    CurrPos_center = CurrPos_multi(:,:,i);
	    CurrPos_center_1 = CurrPos_center(:,1);
	    CurrPos_center_2 = CurrPos_center(:,2);

	    % placeholder for total contact force on the nodes of i-th body
	    cf_contrib_i = zeros(total_nodes, 2);

	    % % Caution: without the trailing transpose ('), for loop misunderstands the sequence as a vector, i.e. without the trailing ', j = [2;3]!!
	    % Must be a row vector for the for loop to sample from.
	    other_particles = (nonzeros( (1:total_particles).*(1:total_particles ~= i) ))';

	    for j = other_particles
		% information about the neighboring particle
		CurrPos_neighbor = CurrPos_multi(:,:,j);
		CurrPos_neighbor_1 = CurrPos_neighbor(:,1);
		CurrPos_neighbor_2 = CurrPos_neighbor(:,2);


		Vol_neighbor = Vol_multi(:,:,j);

		% nodes in the neighboring body that are contact_radius-close to the nodes in the central body
		%% edit this function to output the difference between the points as well!!
		%% distance norm
		%[contact_NbdArr] = gen_NbdArr(CurrPos_center, CurrPos_neighbor, contact_radius, 0);
		[contact_NbdArr] = gen_NbdArr_varlength(CurrPos_center, CurrPos_neighbor, contact_radius, 0);
		% mask for contact neighbors
		mask = ~~contact_NbdArr;

		% volume of the contact neighbors
		cnbd_Vol_neighbor = Vol_neighbor(contact_NbdArr + ~contact_NbdArr) .* mask;

		% direction vectors from the nodes in the central body to the delta-close nodes in the neighboring body
		direction_1 = (CurrPos_neighbor_1(contact_NbdArr + ~contact_NbdArr) - CurrPos_center_1) .* mask;
		direction_2 = (CurrPos_neighbor_2(contact_NbdArr + ~contact_NbdArr) - CurrPos_center_2) .* mask;

		direction_norm = sqrt(direction_1.^2  + direction_2.^2);

		% in the negative direction, i.e. from the neighbors to the center
		direction_unit_1 = - direction_1 ./ (direction_norm + ~direction_norm) .* mask;
		direction_unit_2 = - direction_2 ./ (direction_norm + ~direction_norm) .* mask;

		% flipping the signs
		%% debug
		%direction_unit_1 =  direction_1 ./ (direction_norm + ~direction_norm) .* mask;
		%direction_unit_2 =  direction_2 ./ (direction_norm + ~direction_norm) .* mask;

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


	    %% when contact occurs
	    %if length(contact_NbdArr(16,:)) > 0
	    %i	%=2 is the center
	    %j	%=1 is the neighbor
	    % % contact_NbdArr contains indices of the neighbor 1
	    % % node 16 on the central body, neighbors in neighbor body
	    %contact_NbdArr(16,:)
	    %
	    %disp 'direction_units'
	    %[direction_unit_1(16, :); direction_unit_2(16, :)]
	    %
	    %disp 'dir norm'
	    %direction_norm(16,:)
	    %
	    %disp 'positive contact_radius_contrib'
	    %positive_cont_rad_contrib(16,:)
	    %
	    %fprintf('peridynmic force: (%f, %f)\n', peridynamic_force_multi(16,1,i), peridynamic_force_multi(16,2,i))
	    %
	    %fprintf('contact force: (%f, %f)\n', contact_force_multi(16, 1, i), contact_force_multi(16, 2, i))
	    %fprintf('total contact force: (%f, %f)\n', total_contact_force(16, 1), total_contact_force(16, 2))
	    %
	    %
	    %Quantity = peridynamic_force_multi(:,2,:);
	    %switch save_plot
	    %case 1
	    %savenewpos2_multi(total_particles, CurrPos_multi, Quantity, imgcounter, f, 'multi_', contact_radius)
	    %otherwise
	    %
	    %end
	    %imgcounter = imgcounter + 1;
	    %
	    %pause
	    %end

	end	%endif

    end %endfor loop in i

    % add all the force densities
    total_force_multi = peridynamic_force_multi + contact_force_multi;

    % accelaration
    u0dotdot_multi = (1/rho) .*( total_force_multi + extforce_multi) ;

    %min(u0dotdot_multi(:,2,2))
    %max(u0dotdot_multi(:,2,2))
    %pause

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


    %if (mod(t, 100) == 0) || ( t == 1)
    if mod(t, modulo) == 0
		fprintf('t: %d\n', t)

		location_center = (min(CurrPos_multi(:,2,2)) + max(CurrPos_multi(:,2,2))) /2 * 1e3;
		vel_min = min(u0dot_multi(:,2,2));
		vel_max = max(u0dot_multi(:,2,2));

		store_location(imgcounter) = location_center;
		store_vel_min(imgcounter) = vel_min;
		store_vel_max(imgcounter) = vel_max;

		time = dt * t * 1e3;	% micro seconds

		%fprintf('(i,j) = (%d, %d)\n', i,j)
		%%fprintf('contact neighbors: %d\n', contact_NbdArr(16,:));
		%contact_NbdArr(16,:)

		%fprintf('contact force contribution: %f\n', contact_force_multi(16,:,2))
		%fprintf ('peridynamic force contribution: %f\n', peridynamic_force_multi(16,:,2) )

		% plot this quantity
		%Quantity = contact_force_multi(:,2,:);
		Quantity = peridynamic_force_multi(:,2,:);
		switch save_plot
		case 1
			savenewpos2_multi(total_particles, CurrPos_multi, Quantity, imgcounter, f, 'elastic_collision_', contact_radius, time);
		    otherwise

		end
		imgcounter = imgcounter + 1;
    end


end

toc

% output
NbdArr_out_multi = NbdArr_multi;




