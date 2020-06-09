% generate the mesh

%script

clear all
close all

experiment = 'single_particle' 
%experiment = 'multi_particle'

% plot reference config
plot_reference = 1;
plot_circles = 1;
timesteps = 1800;

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

disp 'Press key to continue'
pause

% generate the neighbor-list
%[NbdArr] = gen_nbdlist(Pos, delta);
%tic
%[NbdArr] = gen_nbdlist2(Pos, delta);	% twice faster
%toc

tic
[NbdArr] = gen_NbdArr(Pos, Pos, delta, 1);
toc

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


disp 'Press key to continue'
pause

% get the position and relative distances of the neighbors
[Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm] = precomputation(NbdArr, Pos, Vol);
[total_nodes, dimension] = size(Pos);

% get the material properties
[delta, rho, Gnot, E, nu, snot, cnot] = material_properties(delta);

% gravity, amplified
gravity = zeros(size(Pos));
gravity(:,2) = - 9.8 * rho * 1e5;

disp 'Precomputation done.'


disp 'Press key to continue'
pause

switch experiment
case 'single_particle'
    % need initial condition
    %[NbdArr_out, u0] = simulate(Pos, NbdArr, Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm);

    % initial data
    uold = zeros(total_nodes,2);
    uolddot = zeros(total_nodes,2);
    uolddotdot = zeros(total_nodes,2);

    % specify
    uolddot = zeros(total_nodes,2) + [1e4 0];

    % external force
    %[extforce] = external_force(geometry, Pos, NbdArr, delta);
    %extforce = gravity;
    extforce = zeros(total_nodes, 2);

    % simulate with initial data
    [NbdArr_out, u0] = simulate_initial(uold, uolddot, uolddotdot, Pos, NbdArr, Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, timesteps);


case 'multi_particle'
    % scaling, shifting, rotation
    %total_particles = 1;
    total_particles = 2;

    % particle location, scaling, and rotation
    particle_scaling = ones(total_particles, 1);
    particle_rotation = zeros(total_particles, 1);
    % % Specify
    particle_shift = [0, 0; 0, 2.5];
    %particle_shift = [0, 0];

    % locations of the particles
    Pos_multi = zeros( [size(Pos), total_particles]);
    for i = 1:total_particles
	Pos_multi(:,:,i) = particle_scaling(i) * Pos + particle_shift(i, :);
    end

    % does volume of individual particle have to change?
    Vol_multi = zeros( [size(Vol), total_particles]);
    for i = 1:total_particles
%% Caution: verify that the volume gets multiplied by the scaling to the power the dimension
	Vol_multi(:,i) = ( (particle_scaling(i)).^dimension) .* Vol;
    end

    % external forces on particles
    extforce_multi = zeros( [size(Pos), total_particles]);
    % gravity applies on the second particle only
    %exforce(:,:,2) = gravity * 1e5;

    % precomputation
    for i = 1:total_particles
	% get the position and relative distances of the neighbors
	[xi_1_multi(:,:,i), xi_2_multi(:,:,i), xi_norm_multi(:,:,i)] = precomputation(NbdArr, Pos, Vol);
    end


    % initial data
     uold_multi = zeros(total_nodes,2, total_particles);
     uolddot_multi = zeros(total_nodes,2, total_particles);
     uolddotdot_multi = zeros(total_nodes,2, total_particles);

     % specify
    uolddot_multi(:,:,2) = zeros(total_nodes,2) + [0 -1e3];
    
    [NbdArr_out, u0_multi] = simulateMultiple(total_particles, uold_multi, uolddot_multi, uolddotdot_multi, Pos_multi, NbdArr, Vol_multi, extforce_multi, delta, rho, cnot, snot, xi_1_multi, xi_2_multi, xi_norm_multi, timesteps);

otherwise
    disp 'No experiment provided'

end

