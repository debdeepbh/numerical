% generate the mesh

%script

clear all
close all


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

scatter(Pos(:,1), Pos(:,2), 5, 'filled')

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
[extforce] = external_force(Pos, NbdArr);
disp 'Press key to continue'
pause


% need initial condition
[NbdArr_out, u0] = simulate2(Pos, NbdArr, Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm);


