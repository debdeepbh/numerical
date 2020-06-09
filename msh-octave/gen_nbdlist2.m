function [NbdArr] = gen_nbdlist2(Pos, delta) 
% Generates the neighborhood array (n = #nodes)
% Input:
%	Pos: Position coordinates of nodes, nx2
%	delta: radius of influence, scalar
%
% Output:
%	NbdArr: nxl matrix, where l is the maximum number of neighbors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Nbd array                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw the peridynamic horizon around the first node
plot_circles = 0

total_nodes = length(Pos);
indices_all = (1:total_nodes)';

% from i-th repetition Pos, subtract i-th row of Pos
diffMat = repmat(Pos, 1, total_nodes) - reshape(Pos', 1, []);

% indices that are within delta distance
% caution: this is a column matrix, i.e. each column has the neighbor indices
neighbor_indices = ( sqrt(diffMat(:, 1:2:end).^2 + diffMat(:, 2:2:end).^2) < delta ) .* indices_all;
% remove the center
neighbor_indices = neighbor_indices .* (neighbor_indices ~= (1:total_nodes));

% max number of neighbors
nbd_count = sum(~~neighbor_indices);
NbdArr = zeros(total_nodes, max(nbd_count));

% remove excess zeros
%% Can we replace this for loop?
for i = 1:total_nodes
	NbdArr(i, 1:nbd_count(i)) = nonzeros( neighbor_indices(:,i)' );
end

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


