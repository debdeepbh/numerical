% generate the mesh

%script

clear all
close all

%experiment = 'single_particle' 
experiment = 'multi_particle'

%material = 'sodalime'
material = 'peridem'

% simulate sodalime crack
simulate_sodalime_prenotch = 0
%simulate_sodalime_prenotch = 1

do_pause = 'no'
%do_pause = 'yes'

% plot reference config
plot_reference = 1;
plot_circles = 1;
%timesteps = 10800;
timesteps = 1e5;	% peridem

modulo = 100;

%geometry = 'sodalime'
%geometry = 'circle'

%delta = 0.002;	% for glass slab
%delta = 0.012;	% peridynamic horizon
%delta = 0.25;	% for the unit circle 
load('delta.mat')

%meshsize = delta/3;
load('meshsize.mat')


disp ('Loading matrices: Pos.mat, Vol.mat, T.mat\n')
load('Pos.mat')
load('Vol.mat')
load('T.mat')
load('geometry.mat')




scatter(Pos(:,1), Pos(:,2), 5, 'filled')
axis equal

switch do_pause
    case 'yes'
	disp 'Press key to continue'
	pause
    otherwise
end

% generate the neighbor-list
%[NbdArr] = gen_nbdlist(Pos, delta);
%tic
%[NbdArr] = gen_nbdlist2(Pos, delta);	% twice faster
%toc

tic
[NbdArr] = gen_NbdArr(Pos, Pos, delta, 1);
toc

if simulate_sodalime_prenotch == 1
    %% for sodalime crack test, borrowing from old code
    load('PosCrack.mat')
    Pos  = PosCrack;
    delta =  0.002;
    dx = delta/4;
    dy = delta/4;
    %Vol = zeros(size(PosCrack)) + dx * dy;
    Vol = zeros(length(PosCrack), 1) + dx * dy;
    %% for sodalime crack test, borrowing from old code
    load('Nbd.mat')
    NbdArr = Nbd;
end

disp('Max neighbors') 
size(NbdArr)

if plot_circles == 1
 hold on
 theta_list = linspace(0, 2* pi , 20);
 %for i = 30:30	%total_nodes
 for i = 1:1	%total_nodes
 	x = Pos(i,1) +  delta * cos(theta_list);
 	y = Pos(i,2) +  delta * sin(theta_list);
 
 	plot(x, y, 'r')

	axis equal
 end
end

switch do_pause
    case 'yes'
	disp 'Press key to continue'
	pause
    otherwise
end


% get the position and relative distances of the neighbors
[Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, nbd_Vol] = precomputation(NbdArr, Pos, Vol);
[total_nodes, dimension] = size(Pos);

% get the material properties
[delta, rho, Gnot, E, nu, snot, cnot, bulk_modulus] = material_properties(delta, material);

% gravity
gravity = zeros(size(Pos));
gravity(:,2) = - 9.8;

disp 'Precomputation done.'

switch do_pause
    case 'yes'
	disp 'Press key to continue'
	pause
    otherwise
end


switch experiment
case 'single_particle'

    % initial data
    uold = zeros(total_nodes,2);
    uolddot = zeros(total_nodes,2);
    uolddotdot = zeros(total_nodes,2);

    extforce = zeros(total_nodes, 2);

%% specify initial data
    %uold = zeros(total_nodes,2) + [0 1e-3];
    uolddot = zeros(total_nodes,2) + [0, 4e4];
    uolddotdot = (zeros(total_nodes, 2) + [0 -16e7] );

    % external force
    %[extforce] = external_force(geometry, Pos, NbdArr, delta);
    extforce = (zeros(total_nodes, 2) + [0 -16e7] ) * rho;

    switch do_pause
	case 'yes'
	    disp 'Press key to continue'
	    pause
	otherwise
    end

    break_bonds = 1;

%dt = 25e-9;
dt = 25e-8;

    % simulate with initial data
    [NbdArr_out, u0] = simulate_initial(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, nbd_Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, timesteps, break_bonds, material, dt);
    %[NbdArr_out, u0] = test_int(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, nbd_Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, timesteps, break_bonds, material, dt);


case 'multi_particle'

    example = 'two_balls'

    switch example
        case 'two_balls'
    	
        otherwise
    	
    end

    % scaling, shifting, rotation
    %total_particles = 1;
    total_particles = 2;


    %contact_radius = delta/2;
    contact_radius = 1.74e-04;	% peridem

    % particle location, scaling, and rotation
    particle_scaling = ones(total_particles, 1);
    particle_rotation = zeros(total_particles, 1);
    % % Specify
    %falling_from = 2.5e-3;	% 5 mm
    falling_from = 5e-3;	% 5 mm
    %virtual_distance = 1e-3 + contact_radius + 1e-3;
    %starting_distance = 2.2e-3;
    starting_distance = 5e-3;
    particle_shift = [0, 0; 0, falling_from];
% particle on the bottom is flipped to ensure symmetric collision
    particle_rotation = [pi; 0];
    %particle_shift = [0, 0];

    % locations of the particles
    Pos_multi = zeros( [size(Pos), total_particles]);

    for i = 1:total_particles
	rot_matrix = [cos(particle_rotation(i)), -sin(particle_rotation(i)); sin(particle_rotation(i)), cos(particle_rotation(i))]
	Pos_multi(:,:,i) = particle_scaling(i) * (rot_matrix * (Pos'))' + particle_shift(i, :);
	%Pos_multi(:,:,i) = particle_scaling(i) * Pos + particle_shift(i, :);
    end

%% rotation
% % rotate the second particle by an angle


    % volume of the nodes
    Vol_multi = zeros( [size(Vol), total_particles]);
    for i = 1:total_particles
%% Caution: verify that the volume gets multiplied by the scaling to the power the dimension
	Vol_multi(:,i) = ( (particle_scaling(i)).^dimension) .* Vol;
    end


    % neighborhood  array
    NbdArr_multi = zeros( [size(NbdArr), total_particles]);
    for i = 1:total_particles
	NbdArr_multi(:,:,i) = NbdArr;
    end

    % precomputation
    for i = 1:total_particles
	% get the position and relative distances of the neighbors
	[xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), nbd_Vol_multi(:,:,i)] = precomputation(NbdArr_multi(:,:,i), Pos_multi(:,:,i), Vol_multi(:,:,i));
    end

    % default initial data
    uold_multi = zeros(total_nodes,2, total_particles);
    uolddot_multi = zeros(total_nodes,2, total_particles);
    uolddotdot_multi = zeros(total_nodes,2, total_particles);

    % default external force density on particles
    extforce_multi = zeros( [size(Pos), total_particles]);

%% specify initial data
    uold_multi(:,:,2) = zeros(total_nodes,2) + [0 (starting_distance - falling_from)];
    uolddot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -60* sqrt(2* 10 * (0.3e-3))];
    %uolddot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -sqrt(2* 10 * (falling_from - starting_distance))];
    %uolddotdot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -10];

  %uolddot_multi(:,:,2) = zeros(total_nodes, 2) + [0, -1.3e-01];	% peridem
  %uolddotdot_multi(:,:,2) = zeros(total_nodes, 2) + [0, -10];	% peridem

    % specify external force density on particles i.e. Force/Volume, or Acceleration * density
    % % gravity applies on the second particle only
    %extforce_multi(:,:,2) = zeros(total_nodes, 2) +  [0, -10 .* rho];


    % % normal stiffness (From Foster's paper)
    %normal_stiffness = 18 * bulk_modulus /( pi * delta^5);
    normal_stiffness = 18 * bulk_modulus /( pi * delta^4);	% from silling-askari05

    % % peridem
    %mass_particle = pi * (1e-3)^2 * rho;
    %meq = mass_particle;% harmonic mean of particles of same size
    %Vmax_by_delmax = (7.385158e+05);
    %normal_stiffness = meq * Vmax_by_delmax^2;

      %normal_stiffness = (7.385158e+05)^2;	% peridem
%% debugging
	%normal_stiffness = normal_stiffness * 1e-5;
    
%dt = 0.02/timesteps;	% peridem
dt = 0.02/1e5;	% peridem

 [NbdArr_out, u0_multi, store_location, store_vel_min, store_vel_max] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr_multi, Vol_multi, nbd_Vol_multi, extforce_multi, normal_stiffness, contact_radius, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, dt, timesteps, delta, modulo);

 time_ss = (1:length(store_location)) *  dt * modulo * 1e3;
 figure
 plot(time_ss, store_location);
 title('Location of the center')
 axis equal

 figure
 plot(time_ss, store_vel_min);
 hold on
 plot(time_ss, store_vel_max);
 hold off
 legend('min', 'max')
 title('Vertical velocity')
 axis equal



otherwise
    disp 'No experiment provided'

end

