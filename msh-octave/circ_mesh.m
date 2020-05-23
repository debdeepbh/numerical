% Define a polygon using

clear all
close all

pkg load msh
pkg load geometry
pkg load fpl

angles = linspace(0, 2*pi, 20)';

P = [ cos(angles), sin(angles)]; 

%  Plot the polygon

drawPolygon(P, '-o')
axis equal

% Generate the `.geo` file from the polygon using `data2geo()` function from `geometry` package:

% this is roughly the mean length of the line segment in the polygon
%meshsize = sqrt (mean (sumsq (diff (P, 1, 1), 2)))/2;
meshsize = 1;
data2geo (P, meshsize, "output", "tempfile.geo");

%Some kind of suggestion to take the mesh size is given by (what is that?)

% Generate the mesh using (without the extension `.geo`. Necessary?)
T = msh2m_gmsh("tempfile");

% plot the mesh
pdemesh (T.p, T.e, T.t);
view(2);
