function [NbdArr] = gen_nbdlist(Pos, delta) 
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


plot_circles = 1

total_nodes = length(Pos);
indices_all = (1:total_nodes)';

nbd_count = zeros(total_nodes,1);
for i = 1:total_nodes
	%neighbor_indices = (sumsq(Pos - Pos(i,:), 2) < delta^2) .* indices_all;
	neighbor_indices = ( sqrt(sumsq(Pos - Pos(i,:), 2))  < delta) .* indices_all;
	% remove the center
	neighbor_indices = neighbor_indices .* (~ (neighbor_indices == i));
	nbdlist{i} = nonzeros(neighbor_indices);

	nbd_count(i) = length(nbdlist{i});
end

% convert to matrix from cell array
NbdArr = zeros(total_nodes, max(nbd_count));
for i = 1:total_nodes
	NbdArr(i, 1:nbd_count(i)) = nbdlist{i};
end

if plot_circles == 1
 hold on
 theta_list = linspace(0, 2* pi , 20);
 for i = 22:22	%total_nodes
 	x = Pos(i,1) +  delta * cos(theta_list);
 	y = Pos(i,2) +  delta * sin(theta_list);
 
 	plot(x, y)
 end
end




