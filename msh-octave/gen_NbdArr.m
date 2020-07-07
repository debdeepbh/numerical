function [NbdArr] = gen_NbdArr(Pos_center, Pos_neighbor, delta, remove_center) 
% Returns the indices of nodes that are within a neighboring body and are within delta distance away from each nodes in the central body 
% Input:
%	Pos_center: Position vectors of the nodes of central body, n x 2
%	Pos_neighbor: Position vectors of the nodes of neighboring body, m x 2
%	delta: specified distance, scalar
%	remove_center: Whether to remove the nodes from which the distance is computed, 0 or 1. Useful when the central body and neighboring body are the same, i.e., generating neighborhood array within the body.
%
% Output:
%	NbdArr: nxl matrix. The i-th row contains the indices of the nodes in the neighboring body which are within delta distance away from the i-th node of the central body. l = max number of neighbors.


[total_nodes, dimension] = size(Pos_center);

indices_all = (1:total_nodes)';

% from i-th repetition of Pos_neighbor, subtract i-th row of Pos
diffMat = repmat(Pos_neighbor, 1, total_nodes) - reshape(Pos_center', 1, []);

% indices that are within delta distance
% caution: this is a column matrix, i.e. each column has the neighbor indices
neighbor_indices = ( sqrt(diffMat(:, 1:dimension:end).^2 + diffMat(:, 2:dimension:end).^2) <= delta ) .* indices_all;

switch remove_center
    case 1
	% remove the center
	neighbor_indices = neighbor_indices .* (neighbor_indices ~= (1:total_nodes));
    otherwise
	
end

% max number of neighbors
nbd_count = sum(~~neighbor_indices);
NbdArr = zeros(total_nodes, max(nbd_count));

% remove excess zeros
%% Can we vectorize this for loop?
for i = 1:total_nodes
	NbdArr(i, 1:nbd_count(i)) = nonzeros( neighbor_indices(:,i)' );
end
