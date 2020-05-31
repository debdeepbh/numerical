% generate the mesh

%script

clear all
close all

single_particle = 1
%multiple_particle = 1

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
[NbdArr] = gen_nbdlist(Pos, delta);

disp('Max neighbors') 
size(NbdArr)
disp 'Press key to continue'
pause

% get the position and relative distances of the neighbors
[Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, NbdVol] = precomputation(NbdArr, Pos, Vol);
disp 'Precomputation done.'

% get the constant external force for the nodes
[extforce] = external_force(geometry, Pos, NbdArr, delta);
disp 'Press key to continue'
pause

if single_particle == 1
    % need initial condition
    [NbdArr_out, u0] = simulate(Pos, NbdArr, Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm);
elseif	multiple_particle == 1

    % simulate multiple particles

    % scaling, shifting, rotation
    total_particles = 3;

    particle_scaling = ones(total_particles, 1);
    particle_rotation = zeros(total_particles, 1);
    particle_shift = [0, 0; 0, 3; 2, 1];

    %[NbdArr_out, u0] = simulateMultiple(Pos, NbdArr, Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm);
end


