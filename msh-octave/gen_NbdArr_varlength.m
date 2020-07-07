function [NbdArr] = gen_NbdArr_varlength(Pos_center, Pos_neighbor, delta, remove_center) 
% Returns the indices of nodes that are within a neighboring body and are within delta distance away from each nodes in the central body 
% Central body and Nbd body can have different number of nodes
% Input:
%	Pos_center: Position vectors of the nodes of central body, n x 2
%	Pos_neighbor: Position vectors of the nodes of neighboring body, m x 2
%	delta: specified distance, scalar
%	remove_center: Whether to remove the nodes from which the distance is computed, 0 or 1. Useful when the central body and neighboring body are the same, i.e., generating neighborhood array within the body.
%
% Output:
%	NbdArr: nxl matrix. The i-th row contains the indices of the nodes in the neighboring body which are within delta distance away from the i-th node of the central body. l = max number of neighbors <= m


[total_nodes_center, dimension] = size(Pos_center);
[total_nodes_neighbor, dimension] = size(Pos_neighbor);

% (i,j)-th element stores the distance^2 between i-th node in the central body and the j-th node in the Neighboring body
% A_ij = |N_j - C_i|
diff = zeros(total_nodes_center, total_nodes_neighbor, dimension);
for i=1:dimension
    % x' - x
    diff(:,:,i) = Pos_neighbor(:,i)' -  repmat(Pos_center(:,i), 1, total_nodes_neighbor) ;
end
% sum in the third dimension to get the norm^2 of the difference
diff_norm = sqrt(sum(diff.^2, 3));
mask = (diff_norm <=delta);
diff_norm_masked =  diff_norm .* mask;

% full list of neighbors, masked
neighbor_indices = repmat(1:total_nodes_neighbor, total_nodes_center, 1) .* mask;

switch remove_center
    case 1
	% remove the center
	neighbor_indices = neighbor_indices .* (neighbor_indices ~= (1:total_nodes_center)' );
    otherwise
	
end

% remove excess zeros
nbd_count = sum(~~neighbor_indices, 2);
max_nbd = max(nbd_count);
NbdArr = zeros(total_nodes_center, max_nbd);

%% Can we vectorize this for loop?
for i = 1:total_nodes_center
	NbdArr(i, 1:nbd_count(i)) = nonzeros( neighbor_indices(i,:) );
end
