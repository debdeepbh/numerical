% generate the mesh
script

% column vec
Pos = T_centroid';
delta = 0.7;

NbdArr = gen_nbdlist(Pos, delta);
