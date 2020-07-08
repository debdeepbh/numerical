function [NbdArr_out_multi, u0_multi] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr_multi, Vol_multi, nbd_Vol_multi, extforce_multi, normal_stiffness, contact_radius, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, dt, timesteps, delta, modulo, break_bonds, with_wall, allow_friction, friction_coefficient);


% break bonds or not

%save_plot = 0;
save_plot = 1;

use_influence_function = 0
%use_influence_function = 1

close all

[total_nodes, temp, temp] = size(NbdArr_multi);

switch with_wall
    case 1
	load('wall_Pos');
	load('wall_Vol');
	load('wall_T');

	% wall dimension
	load('wall_left');
	load('wall_right');
	load('wall_top');
	load('wall_bottom');

	% wall displacement
	wall_u0 = zeros(size(wall_Pos));
    otherwise
	
end

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

    switch with_wall
        case 1
	    % current wall position
	    wall_CurrPos = wall_u0 + wall_Pos ;
        otherwise
    end


    % placeholder for contact forces on nodes due to neighboring particles
    peridynamic_force_multi = zeros(size(u0_multi));
    contact_force_multi = zeros(size(u0_multi));
    friction_force_multi = zeros(size(u0_multi));

%% for debugging, only one particle
    %for i=2
    for i=1:total_particles

	% information about the central particle
	CurrPos_center = CurrPos_multi(:,:,i);
	min_Pos_center = min(CurrPos_center);
	max_Pos_center = max(CurrPos_center);
	Vol_center = Vol_multi(:,i);

	vel_center_1 = uolddot_multi(:,1,i);
	vel_center_2 = uolddot_multi(:,2,i);

	% internal force - peridynamic
	[peridynamic_force_multi(:,:,i), stretch_multi(:,:,i)] = peridynamic_force_bypos(CurrPos_center, NbdArr_multi(:,:,i), nbd_Vol_multi(:,:,i), xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), cnot, delta, use_influence_function) ;

	% pairwise contact forces
%% this if statement can probably be safely removed, since other_particles will produce empty list when total_particles == 1
	if total_particles > 1
	    % % Caution: without the trailing transpose ('), for loop misunderstands the sequence as a vector, i.e. without the trailing ', j = [2;3]!!
	    % Must be a row vector for the for loop to sample from.
	    other_particles = (nonzeros( (1:total_particles).*(1:total_particles ~= i) ))';

	    % placeholder for various types of forces on the nodes of i-th body
	    cf_contrib_on_center = zeros(total_nodes, 2);
	    ff_contrib_on_center = zeros(total_nodes, 2);
	    for j = other_particles
		% information about the neighboring particle
		CurrPos_neighbor = CurrPos_multi(:,:,j);
		min_Pos_neighbor = min(CurrPos_neighbor);
		max_Pos_neighbor = max(CurrPos_neighbor);

		if intersects_box(min_Pos_center - contact_radius, max_Pos_center + contact_radius, min_Pos_neighbor, max_Pos_neighbor)

		    fprintf('contact (%d, %d) \n', i,j);

		    Vol_neighbor = Vol_multi(:,j);
		    % nodes in the neighboring body that are contact_radius-close to the nodes in the central body
		    %% edit this function to output the difference between the points as well!!
		    [contact_NbdArr] = gen_NbdArr_varlength(CurrPos_center, CurrPos_neighbor, contact_radius, 0);

		    % get the contact for on the center due to the neighbor
		    [total_contact_force, nbd_contact_force_1, nbd_contact_force_2, nbd_direction_unit_1, nbd_direction_unit_2] = contact_force(contact_NbdArr, CurrPos_center, CurrPos_neighbor, Vol_neighbor, contact_radius, normal_stiffness);

		    %% add this force to total internal force
		    cf_contrib_on_center = cf_contrib_on_center + total_contact_force;

		    % % Tangential friction force
		    total_friction_force = zeros(size(Pos_multi));
		    switch allow_friction
		    case 1
			% contact velocity
			% % Caution: this should be computed from u0dot actually, but that is not available until the end of the loop
			vel_neighbor_1 = uolddot_multi(:,1,j);
			vel_neighbor_2 = uolddot_multi(:,2,j);

			% compute the total friction force acting on the nodes in the central body due the neighboring body
			[total_friction_force] = friction_force(contact_NbdArr, vel_center_1, vel_center_2, vel_neighbor_1, vel_neighbor_2, nbd_contact_force_1, nbd_contact_force_2, nbd_direction_unit_1, nbd_direction_unit_2, Vol_center, friction_coefficient);

			% add the friction force contribution from neighbor j
			ff_contrib_on_center = ff_contrib_on_center + total_friction_force;
		    otherwise
		    end

		end	%endif close-by

	    end	% endfor loop in j

	    contact_force_multi(:,:, i) = cf_contrib_on_center;
	    friction_force_multi(:,:, i) = ff_contrib_on_center;
	end	%endif total_particles > 1

	% % interaction between the i-th particle (center) and the wall
	switch with_wall
	case 1
	    % if near the wall
	    if (min_Pos_center(1) < (wall_left + contact_radius)) | ( min_Pos_center(2) < (wall_bottom + contact_radius) ) | (max_Pos_center(1) > (wall_right - contact_radius)) | (max_Pos_center(2) > (wall_top - contact_radius))

		%fprintf('wall contact %d\n',i);

		[wall_contact_NbdArr] = gen_NbdArr_varlength(CurrPos_center, wall_CurrPos, contact_radius, 0);
		%[wall_contact_force] = contact_force(wall_contact_NbdArr, CurrPos_center, wall_CurrPos, wall_Vol, contact_radius, normal_stiffness) ;
		[wall_contact_force, wall_nbd_contact_force_1, wall_nbd_contact_force_2, wall_nbd_direction_unit_1, wall_nbd_direction_unit_2] = contact_force(wall_contact_NbdArr, CurrPos_center, wall_CurrPos, wall_Vol, contact_radius, normal_stiffness);

		% add it to the total contact force on the i-th particle due to wall
		contact_force_multi(:,:, i) = contact_force_multi(:,:,i) + wall_contact_force;
		switch allow_friction
		case 1
			% contact velocity

			% velocity of neighbor, i.e., the wall
			vel_neighbor_1 = zeros(size(wall_Pos));
			vel_neighbor_2 = zeros(size(wall_Pos));

			% compute the total friction force acting on the nodes in the central body due the neighboring body
			[wall_friction_force] = friction_force(wall_contact_NbdArr, vel_center_1, vel_center_2, vel_neighbor_1, vel_neighbor_2, wall_nbd_contact_force_1, wall_nbd_contact_force_2, wall_nbd_direction_unit_1, wall_nbd_direction_unit_2, Vol_center, friction_coefficient);

			% add the friction force contribution from neighbor j
			friction_force_multi(:,:,i) = friction_force_multi(:,:,i) + wall_friction_force;
		otherwise
		end



	    end
	otherwise
	end

	%% end of contact for computations

    end %endfor loop in i

    % add all the force densities
    total_force_multi = peridynamic_force_multi + contact_force_multi + friction_force_multi;

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

	file_string = 'no_friction_collision_';

	fprintf('t: %d\n', t)

	time = dt * t * 1e3;	% micro seconds

	% plot this quantity
	%Quantity = contact_force_multi(:,2,:);
	Quantity = peridynamic_force_multi(:,2,:);
	switch save_plot
	case 1
	    switch with_wall
	    case 1
		savenewpos2_multi_with_wall(total_particles, CurrPos_multi, wall_CurrPos, Quantity, imgcounter, f, file_string, contact_radius, time);
	    otherwise
		savenewpos2_multi(total_particles, CurrPos_multi, Quantity, imgcounter, f, file_string, contact_radius, time);
	    end
	otherwise

	end
	imgcounter = imgcounter + 1;
    end


end

toc

% output
NbdArr_out_multi = NbdArr_multi;




