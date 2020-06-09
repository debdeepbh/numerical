function [Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm] = precomputation(NbdArr, Pos, Vol) 
% Computes the Position and relative distance (xi, xi_norm) of the neighbors from the center of the peridynamic circle. n = #nodes, l = max number of neighbors
% Input:
%	NbdArr: matrix containing the indices of the neighbors, nxl
%	Pos: position of the nodes nx2
%
% Output:
%	Nbd_Pos_1: x-position of the neighbors, nxl
%	Nbd_Pos_2: x-position of the neighbors, nxl
%	xi_1: x-value of the distance vector from the center of the peridynamic circle to the neighbor, nxl
%	xi_2: y-value of the distance vector from the center of the peridynamic circle to the neighbor, nxl
%	xi_norm: distance (scalar) of the neighbor from the center of the peridynamic circle, nxl
%	NbdVol: nodal volume associated with each of the neighbors of the nodes, nxl


% masking matrix
Mask = (NbdArr > 0) ;

Pos_1 = Pos(:,1);
Pos_2 = Pos(:,2);

% position of the neighbors
Nbd_Pos_1 = Pos_1(NbdArr + ~NbdArr) .* Mask;
Nbd_Pos_2 = Pos_2(NbdArr + ~NbdArr) .* Mask;

% xi associated to the neighbor elements
Nbd_xi_1 = (Nbd_Pos_1 - Pos_1) .* Mask;
Nbd_xi_2 = (Nbd_Pos_2 - Pos_2) .* Mask;

Nbd_xi_norm = sqrt(Nbd_xi_1.^2 + Nbd_xi_2.^2);

% nodal volume
%NbdVol = Vol(NbdArr + ~NbdArr) .* (~~NbdArr);
