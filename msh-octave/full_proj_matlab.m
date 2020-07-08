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

% break bonds or not
break_bonds = 0;
%break_bonds = 1

%with_wall = 0
with_wall = 1

allow_friction = 0
%allow_friction = 1
friction_coefficient = 0.5;

%total_particles = 1;
%total_particles = 3;	% 3 particles

%specified_initial_data = '3_equidistant'
%specified_initial_data = '2_vertical'
%specified_initial_data = 'friction_test'
specified_initial_data = 'falling_tube'

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
%tic
%[NbdArr] = gen_NbdArr(Pos, Pos, delta, 1);	% fastest, not anymore
%toc
tic
[NbdArr] = gen_NbdArr_varlength(Pos, Pos, delta, 1);	% fastest
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


%dt = 25e-9;
dt = 25e-8;

    % simulate with initial data
    [NbdArr_out, u0] = simulate_initial(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, nbd_Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, timesteps, break_bonds, material, dt);
    %[NbdArr_out, u0] = test_int(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, nbd_Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, timesteps, break_bonds, material, dt);


case 'multi_particle'


    dt = 0.02/1e5;	% peridem

    %contact_radius = delta/2;
    contact_radius = 1.74e-04;	% peridem
    % % normal stiffness (From Foster's paper)
    %normal_stiffness = 18 * bulk_modulus /( pi * delta^5);
    normal_stiffness = 18 * bulk_modulus /( pi * delta^4);	% from silling-askari05

    % % Specify
    switch specified_initial_data
    case 'falling_tube'
	falling_from = 3e-3;	% 5 mm
	%starting_distance = 2.2e-3;
	starting_distance = 3e-3;
	particle_shift = [0, falling_from];	% 2 particles
	%particle_rotation = [pi; 0];	% for elastic collision of unit circles % 2 particles

	%% Initial data
	uold_multi(:,:,1) = zeros(total_nodes,2) + [0 (starting_distance - falling_from)];
	uolddot_multi(:,:,1) = zeros(total_nodes,2) +  [0, -60* sqrt(2* 10 * (0.3e-3))];
	%uolddotdot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -10];
    case 'friction_test'
	falling_from = 3e-3;	% 5 mm
	%starting_distance = 2.2e-3;
	starting_distance = 3e-3;
	particle_shift = [0, 0; -2e-3 + contact_radius, falling_from];	% 2 particles
	%particle_shift = [0, 0; 0, falling_from];	% 2 particles
	%particle_rotation = [pi; 0];	% for elastic collision of unit circles % 2 particles
	particle_rotation = [pi/6; -pi/6];

	%% Initial data
	uold_multi(:,:,2) = zeros(total_nodes,2) + [0 (starting_distance - falling_from)];
	uolddot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -60* sqrt(2* 10 * (0.3e-3))];
	%uolddotdot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -10];
    case '2_vertical'
	%falling_from = 2.5e-3;	% 5 mm
	falling_from = 3e-3;	% 5 mm
	%starting_distance = 2.2e-3;
	starting_distance = 3e-3;
	particle_shift = [0, 0; 0, falling_from];	% 2 particles

	% particle on the bottom is flipped to ensure symmetric collision
	%particle_rotation = [pi; 0];	% for elastic collision of unit circles % 2 particles
	particle_rotation = [pi - pi/8; -pi/8];	% for pacman collision
	%particle_rotation = [pi - pi/8; -pi/8];	% for circle_w_prenotch
	%particle_rotation = [0; 0];	% for triangles

	%% Initial data
	uold_multi(:,:,2) = zeros(total_nodes,2) + [0 (starting_distance - falling_from)];
	uolddot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -60* sqrt(2* 10 * (0.3e-3))];	% 2 particles
	%uolddot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -sqrt(2* 10 * (falling_from - starting_distance))];
	uolddotdot_multi(:,:,2) = zeros(total_nodes,2) +  [0, -10];
	%uolddot_multi(:,:,2) = zeros(total_nodes, 2) + [0, -1.3e-01];	% peridem
	%uolddotdot_multi(:,:,2) = zeros(total_nodes, 2) + [0, -10];	% peridem

	% specify external force density on particles i.e. Force/Volume, or Acceleration * density
	% % gravity applies on the second particle only
	%extforce_multi(:,:,2) = zeros(total_nodes, 2) +  [0, -10 .* rho];

    case '3_equidistant'
	particle_shift = [-5/2, 0; 5/2, 0; 0, 5 * sin(pi/3)] * 1e-3; % 3 particles
	particle_rotation = [pi/3; pi/3; pi/3]; % 3 particles

	%% initial data
	uolddot_multi(:,:,1) = zeros(total_nodes,2) +  -[sqrt(3)/2, 1/2 ] * -60* sqrt(2* 10 * (0.3e-3));	% 3 particles
	uolddot_multi(:,:,2) = zeros(total_nodes,2) +  -[-sqrt(3)/2, 1/2 ] * -60* sqrt(2* 10 * (0.3e-3));	% 3 particles
	uolddot_multi(:,:,3) = zeros(total_nodes,2) +  -[0, -1] * -60* sqrt(2* 10 * (0.3e-3));	% 3 particles

    otherwise
	disp('No initial setup specified, assuming default.')
    end

    %% Necessary variables, if not defined
    %% geometry
    if ~exist('total_particles', 'var')
	% compute the total number of particles from the shift
	[total_particles, temp] = size(particle_shift);
	total_particles
    else
	fprintf('specified total_particles = %d',total_particles);
    end
    if ~exist('particle_scaling', 'var')
	% particle location, scaling, and rotation
	particle_scaling = ones(total_particles, 1);
    end
    if ~exist('particle_rotation', 'var')
	particle_rotation = zeros(total_particles, 1);
    end
    if ~exist('particle_shift', 'var')
	particle_shift = zeros(total_particles, dimension);
    end

    %% initial data
    if ~exist('uold_multi', 'var')
	uold_multi = zeros(total_nodes,2, total_particles);
    end
    if ~exist('uolddot_multi', 'var')
	uolddot_multi = zeros(total_nodes,2, total_particles);
    end
    if ~exist('uolddotdot_multi', 'var')
	uolddotdot_multi = zeros(total_nodes,2, total_particles);
    end
    if ~exist('extforce_multi', 'var')
	% default external force density on particles
	extforce_multi = zeros( [size(Pos), total_particles]);
    end

    % generate the geometry
    %Pos_multi = zeros( [size(Pos), total_particles]);
    %Vol_multi = zeros( [size(Vol), total_particles]);
    %NbdArr_multi = zeros( [size(NbdArr), total_particles]);
    for i = 1:total_particles
	rot_matrix = [cos(particle_rotation(i)), -sin(particle_rotation(i)); sin(particle_rotation(i)), cos(particle_rotation(i))]
	Pos_multi(:,:,i) = particle_scaling(i) * (rot_matrix * (Pos'))' + particle_shift(i, :);
	Vol_multi(:,i) = ( (particle_scaling(i)).^dimension) .* Vol;
	NbdArr_multi(:,:,i) = NbdArr;
	%% precomputation
	% get the position and relative distances of the neighbors
	[xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i), nbd_Vol_multi(:,:,i)] = precomputation(NbdArr_multi(:,:,i), Pos_multi(:,:,i), Vol_multi(:,i));
    end


 [NbdArr_out, u0_multi] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr_multi, Vol_multi, nbd_Vol_multi, extforce_multi, normal_stiffness, contact_radius, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, dt, timesteps, delta, modulo, break_bonds, with_wall, allow_friction, friction_coefficient);

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

