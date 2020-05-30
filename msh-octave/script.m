
clear all
close all
% Define a polygon using

%P = [10 10; 40 15; 40 35; 20 25; 15 30];
P = [0 0; 0.1 0; 0.1 0.04; 0 0.04];

% circle, Caution: the last node should _NOT_ overlap with the first!
angles = linspace(0, 2*pi- 2*pi/10, 10)';
%P = [ cos(angles), sin(angles)]; 

%meshsize = sqrt (mean (sumsq (diff (P, 1, 1), 2)))/2/8;
%meshsize = 0.5;
meshsize = 1/64;

% Load the geom package

pkg load geometry

%  Plot the polygon

drawPolygon(P, '-o')

pause (1)

% Generate the `.geo` file from the polygon using `data2geo()` function from `geometry` package:


data2geo (P, meshsize, "output", "tempfile.geo");

%Some kind of suggestion to take the mesh size is given by (what is that?)


% Load the meshing library
pkg load msh
% Generate the mesh using (without the extension `.geo`. Necessary?)
T = msh2m_gmsh("tempfile");

% Load the FEM plotting library `fpl`
pkg load fpl

% plot the mesh
pdemesh (T.p, T.e, T.t);
view(2);


% plot the node indices
text(T.p(1,:), T.p(2,:), num2cell(1:length(T.p), 1) )

% plot the element indices
[T_centroid, T_area] = msh2m_geometrical_properties(T, "bar", "area");
text(T_centroid(1,:), T_centroid(2,:), num2cell(1:length(T_centroid), 1) )

hold on
scatter(T_centroid(1,:), T_centroid(2,:), T_area'/min(T_area)*300)

pause (1)

%%% Topological properties
% Get the element indices of neighbors of each element (`NaN` for missing values)
% msh2m_topological_properties(T, "n")

% bounday elements
T_boundary = msh2m_topological_properties(T, "boundary");

axis equal

