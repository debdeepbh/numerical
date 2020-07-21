function [NbdArr_out_multi, u0_multi] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr_multi, Vol_multi, nbd_Vol_multi, extforce_multi, normal_stiffness, contact_radius, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, dt, timesteps, delta, modulo, break_bonds, with_wall, allow_friction, friction_coefficient, allow_contact, allow_damping, damping_ratio, wall_type);

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

	% wall initial data
	wall_u0 = zeros(size(wall_Pos));
	wall_u0dot = zeros(size(wall_Pos));
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
    damping_force_multi = zeros(size(u0_multi));

%% for debugging, only one particle
    %for i=2
    for i=1:total_particles

	% information about the central particle
	CurrPos_center = CurrPos_multi(:,:,i);
	min_Pos_center = min(CurrPos_center);
	max_Pos_center = max(CurrPos_center);
	Vol_center = Vol_multi(:,i);
	velocity_center = uolddot_multi(:,:,i);

	% internal force - peridynamic
	[peridynamic_force_multi(:,:,i), stretch_multi(:,:,i)] = peridynamic_force_bypos(CurrPos_center, NbdArr_multi(:,:,i), nbd_Vol_multi(:,:,i), xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), cnot, delta, use_influence_function) ;

	% % pairwise forces

	% placeholder for various types of forces on the nodes of i-th body
	cf_contrib_on_center = zeros(total_nodes, 2);
	ff_contrib_on_center = zeros(total_nodes, 2);
	df_contrib_on_center = zeros(total_nodes, 2);

	% % Caution: without the trailing transpose ('), for loop misunderstands the sequence as a vector, i.e. without the trailing ', j = [2;3]!!
	% Must be a row vector for the for loop to sample from.
	other_particles = (nonzeros( (1:total_particles).*(1:total_particles ~= i) ))';
	for j = other_particles
	    % information about the neighboring particle
	    CurrPos_neighbor = CurrPos_multi(:,:,j);
	    min_Pos_neighbor = min(CurrPos_neighbor);
	    max_Pos_neighbor = max(CurrPos_neighbor);

	    % if nearby
	    if intersects_box(min_Pos_center - contact_radius, max_Pos_center + contact_radius, min_Pos_neighbor, max_Pos_neighbor)

		%fprintf('contact (%d, %d) \n', i,j);

		Vol_neighbor = Vol_multi(:,j);
		velocity_neighbor = uolddot_multi(:,:,j);

		[total_contact_force, total_friction_force, total_damping_force] = pairwise_forces(CurrPos_center, CurrPos_neighbor, velocity_center, velocity_neighbor, Vol_neighbor, contact_radius, normal_stiffness, friction_coefficient, damping_ratio, rho);

	%disp('damping')
	%max(max(total_contact_force))
	%max(max(total_damping_force))

		cf_contrib_on_center = cf_contrib_on_center + total_contact_force;
		ff_contrib_on_center = ff_contrib_on_center + total_friction_force;
		df_contrib_on_center = df_contrib_on_center + total_damping_force;
	    end	%endif close-by
	end	% endfor loop in j

	contact_force_multi(:,:, i) = cf_contrib_on_center;
	friction_force_multi(:,:, i) = ff_contrib_on_center;
	damping_force_multi(:,:, i) = df_contrib_on_center;

	% % interaction between the i-th particle (center) and the wall
	if with_wall
	    wall_min = [wall_left, wall_bottom];
	    wall_max = [wall_right, wall_top];

	    switch wall_type
	    case 'box'
		% % for interior of a box
		wall_contact_true = ~within_interior(min_Pos_center - contact_radius, max_Pos_center + contact_radius, wall_min, wall_max);
	    case 'rectangle'
		% % for rectangular object
		wall_contact_true =  (min_Pos_center(1) < (wall_left + contact_radius)) | ( min_Pos_center(2) < (wall_bottom + contact_radius) ) | (max_Pos_center(1) > (wall_right - contact_radius)) | (max_Pos_center(2) > (wall_top - contact_radius));
	    otherwise
		disp 'wrong wall type'
		return
	    end	%endswitch wall_type

	    if wall_contact_true 
		%fprintf('wall contact %d\n',i);


		% here, neighbor is the wall
		[wall_contact_force, wall_friction_force, wall_damping_force] = pairwise_forces(CurrPos_center, wall_CurrPos, velocity_center, wall_u0dot, wall_Vol, contact_radius, normal_stiffness, friction_coefficient, damping_ratio, rho);

		friction_force_multi(:,:,i) = friction_force_multi(:,:,i) + wall_friction_force;
		contact_force_multi(:,:,i) = contact_force_multi(:,:,i) + wall_contact_force;
		damping_force_multi(:,:,i) = damping_force_multi(:,:,i) + wall_damping_force;
	    end
	end %endif with_wall

    end %endfor loop in i

    % add all the force densities
    total_force_multi = peridynamic_force_multi;
    if allow_contact == 1
	total_force_multi = total_force_multi + contact_force_multi;
    end
    if allow_friction == 1
	total_force_multi = total_force_multi + friction_force_multi;
    end
    if allow_damping == 1
	total_force_multi = total_force_multi + damping_force_multi;
    end

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

    %if (mod(t, 100) == 0) || ( t == 1)
    if mod(t, modulo) == 0

	file_string = 'many_particles_';

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


%     %% Store quantities
%    store_time(t,:) = dt * t * 1e3;	% micro s
%    store_center(t,:) = mean(CurrPos_multi(:,:,2));
%    store_CurrPos_min(t,:) = min(CurrPos_multi(:,:,2));
%    %store_CurrPos_max(t,:) = max(CurrPos_multi(:,:,2));
    

end

toc

% output
NbdArr_out_multi = NbdArr_multi;



