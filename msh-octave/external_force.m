function [extforce] = external_force(Pos, NbdArr) 
% Specify force (b) on the nodes of the body, n=#nodes
% Input:
%	: 
%
% Output:
%	extforce: Force vector, nx2

% plot force
plot_force = 1;

gravity = 9.8;
totalnodes = length(Pos);
extforce = zeros(totalnodes,2); % extforce(i,:) = external force on i-th serial

% %%%%%%%%%%%%%%5



% get the placeholders of the nodes on the top and bottom edges of the rectangle
top_edge_location = (abs(Pos(:,2) - 0.04) < 1e-10);
bottom_edge_location = (abs(Pos(:,2) - 0) < 1e-10);


% distribute the external force over the 
%sigma = 12e6/sum(top_edge_location);
%dx = delta/4;
dx = 0.002/4;
%sigma = 12e6/sum(top_edge_location);
sigma = 12e6/dx;

% hopefully the top and bottom edges are non-overlapping
extforce = [zeros(totalnodes, 1), top_edge_location] .* sigma + [zeros(totalnodes, 1), bottom_edge_location] .* (-sigma);

if plot_force == 1
	disp 'Plotting the force'
	hold on
	quiver(Pos(:,1), Pos(:,2), extforce(:,1), extforce(:,2))
end