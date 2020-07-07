function [NbdArr] = gen_NbdArr_varlength(Pos_center, Pos_neighbor, delta, remove_center) 
% Central body and Nbd body can have different number of nodes
% Returns the indices of nodes that are within a neighboring body and are within delta distance away from each nodes in the central body 
% Input:
%	Pos_center: Position vectors of the nodes of central body, n x 2
%	Pos_neighbor: Position vectors of the nodes of neighboring body, m x 2
%	delta: specified distance, scalar
%	remove_center: Whether to remove the nodes from which the distance is computed, 0 or 1. Useful when the central body and neighboring body are the same, i.e., generating neighborhood array within the body.
%
% Output:
%	NbdArr: nxl matrix. The i-th row contains the indices of the nodes in the neighboring body which are within delta distance away from the i-th node of the central body. l = max number of neighbors.


[total_nodes_center, dimension] = size(Pos_center);
[total_nodes_neighbor, dimension] = size(Pos_neighbor);

% (i,j)-th element stores the distance^2 between i-th node in the central body and the j-th node in the Neighboring body
% A_ij = |C_i - N_j|
sum_diff_sq = zeros(total_nodes_center, total_nodes_neighbor);
for i=1:dimension
    diff = repmat(Pos_center(:,i), 1, total_nodes_neighbor) - Pos_neighbor(:,i)';
    sum_diff_sq = sum_diff_sq + diff.^2;
end

neighbor_indices = (sqrt(sum_diff_sq) <= delta) .* repmat(1:total_nodes_neighbor, total_nodes_center, 1);


switch remove_center
    case 1
	% remove the center
	neighbor_indices = neighbor_indices .* (neighbor_indices ~= (1:total_nodes_center)' );
    otherwise
	
end


% number of neighbors
nbd_count = sum(~~neighbor_indices, 2);
max_nbd = max(nbd_count);
NbdArr = zeros(total_nodes_center, max_nbd);

% remove excess zeros
%% Can we vectorize this for loop?
for i = 1:total_nodes_center
	NbdArr(i, 1:nbd_count(i)) = nonzeros( neighbor_indices(i,:) );
end
