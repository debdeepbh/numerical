% generate the mesh

%script

clear all
close all


geometry='peridem'
%geometry='circle'
%geometry='sodalime'

switch geometry
case 'sodalime'
	%delta = 0.002;	% for glass slab
	delta = 0.004;	% for glass slab
	%delta = 0.012;	% peridynamic horizon
	%delta = 0.25;	% for the unit circle 
	meshsize = delta/2;
case 'peridem'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = 1e-3/5;
otherwise

end

%[Pos, BdryNodes, Vol] = genmesh(delta);
[Pos, Vol, T] = genmesh2(geometry, meshsize);

printf('total nodes: %d \n', length(Pos));

printf('Saving matrices: Pos.mat, Vol.mat, T.mat\n')

% save in matlab-readable format
save -mat7-binary 'Pos.mat' 'Pos'
save -mat7-binary 'Vol.mat' 'Vol'
save -mat7-binary 'T.mat' 'T'
save -mat7-binary 'delta.mat' 'delta'
save -mat7-binary 'meshsize.mat' 'meshsize'
save -mat7-binary 'geometry.mat' 'geometry'

disp 'Press key to continue'
pause

% generate the neighbor-list
[NbdArr] = gen_nbdlist(Pos, delta);

[nn, l] = size(NbdArr);

printf('Max neighbors: %d \n', l);
disp 'Press key to continue'
pause

% get the position and relative distances of the neighbors
[xi_1, xi_2, xi_norm, NbdVol] = precomputation(NbdArr, Pos, Vol);
disp 'Precomputation done.'

% get the constant external force for the nodes
[extforce] = external_force(Pos, NbdArr);
disp 'Press key to continue'
pause


% need initial condition
[NbdArr_out, u0, u_norm, disp] = simulate(Pos, NbdArr, NbdVol, extforce, delta, xi_1, xi_2, xi_norm);


