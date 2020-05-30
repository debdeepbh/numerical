function [Pos, Vol, T] = genmesh2(meshsize) 
% Generates a mesh with boundary
% Input:
%	meshsize: 
%
% Output:
%	Pos:
%	BdryNodes:
%	Vol:

% plot the indices of the nodes, very slow when there are many nodes
draw_text = 0;

% Define a polygon using
%P = [10 10; 40 15; 40 35; 20 25; 15 30];
P = [0.1 0; 0.2 0; 0.2 0.04;0.1 0.04];

% circle, Caution: the last node should _NOT_ overlap with the first!
steps = 20;
angles = linspace(0, 2*pi- 2*pi/steps, steps)';
%P = [ cos(angles), sin(angles)]; 

% Load the geom package

pkg load geometry

%%  Plot the polygon
%
%drawPolygon(P, '-o')
%pause (1)

% Generate the `.geo` file from the polygon using `data2geo()` function from `geometry` package:


data2geo (P, meshsize, "output", "tempfile.geo");


tic

% Load the meshing library
pkg load msh
% Generate the mesh using (without the extension `.geo`. Necessary?)
T = msh2m_gmsh("tempfile");

toc

% Load the FEM plotting library `fpl`
pkg load fpl


% plot the mesh
pdemesh (T.p, T.e, T.t);
view(2);

% row represents the node, each row contains the indices of the triangles it shares
node_nbd_element = element2node(T.t(1:3,:)');




if draw_text == 1
	% plot the node indices
	text(T.p(1,:), T.p(2,:), num2cell(1:length(T.p), 1) );
end

% get the centroid and area of the mesh elements
[T_centroid, T_area] = msh2m_geometrical_properties(T, "bar", "area");


% distribute the volume of the triangles equally to all 3 nodes
node_area_distributed = T_area(node_nbd_element + ~node_nbd_element) .* (~~node_nbd_element) /3;
% sum up the volume contribution from all neighboring elements
Vol = sum(node_area_distributed, 2);

%%% Topological properties
% Get the element indices of neighbors of each element (`NaN` for missing values)
% msh2m_topological_properties(T, "n")

axis equal

% column vec
Pos = (T.p)';


