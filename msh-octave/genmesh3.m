
% generate the mesh

%script

clear all
close all


%geometry='peridem'
%geometry='pacman'
%geometry='peridem_wall'
geometry='peridem_tube'
%geometry = 'peridem_w_prenotch'
%geometry = 'peridem_triangle'
%geometry='unitcircle'
%geometry='sodalime'


draw_text = 1;

switch geometry
case 'sodalime'
	%delta = 0.002;	% for glass slab
	delta = 0.004;	% for glass slab
	%delta = 0.012;	% peridynamic horizon
	%delta = 0.25;	% for the unit circle 
	meshsize = delta/2;
	% Define a polygon using
	%P = [0.1 0; 0.2 0; 0.2 0.04; 0.1 0.04];

	% with pre-notch
	P = [0.1, 0; 0.2, 0; 0.2, 0.04; 0.1, 0.04; 0.1, (0.04/2) + meshsize/2; (0.1 + 0.2)/2, (0.04/2) + meshsize/2; (0.1 + 0.2)/2, (0.04/2) - meshsize/2; 0.1, (0.04/2) - meshsize/2 ];
	%%  Plot the polygon
	%
	%drawPolygon(P, '-o')
	%pause (1)
case 'peridem'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% circle, Caution: the last node should _NOT_ overlap with the first!
	steps = 20;
	angles = linspace(0, 2*pi- 2*pi/steps, steps)';
	P = 1e-3 * [ cos(angles), sin(angles)]; % 1mm size
case 'pacman'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% circle, Caution: the last node should _NOT_ overlap with the first!
	steps = 20;
	% opening of 90 degrees
	angle_mouth = pi/2;
	angles = linspace(0, 2*pi- angle_mouth, steps)';
	P = 1e-3 * [ cos(angles), sin(angles)]; % 1mm size
	P = [P; 0, 0]
case 'unitcircle'
	meshsize = 1/5;
	delta = meshsize;	% peridem, all the nodes are neighbors
	% circle, Caution: the last node should _NOT_ overlap with the first!
	steps = 20;
	angles = linspace(0, 2*pi- 2*pi/steps, steps)';
	P = [ cos(angles), sin(angles)]; % 1mm size
case 'peridem_w_prenotch'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% circle, Caution: the last node should _NOT_ overlap with the first!
	steps = 20;
	% opening of 90 degrees
	angle_mouth = pi/2;
	angles = linspace(0 + 2*pi/steps/2, 2*pi- 2*pi/steps/2, steps)';
	P = 1e-3 * [ cos(angles), sin(angles)]; % 1mm size
	P = [P; 0, 0]
    case 'peridem_triangle'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;

	side_length = 2e-3;
	P = [-side_length/2, 0; side_length/2 , 0; 0, side_length * sin(pi/3)]
    case 'peridem_wall'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	%a = 5*meshsize * 1e3;	% in mm
	a = meshsize * 1e3;	% in mm
		% for rectangular wall
		wall_bottom = -2 * 1e-3;
		wall_left = -2 * 1e-3;
		wall_right = 2 * 1e-3;
		wall_top = 10;	% fake
	P = 1e-3 * [-2,3; -2,-2; 2,-2; 2,3; 2+a,3; 2+a,-2-a; -2-a, -2-a; -2-a,3];
	%P = 1e-3 * [-2,3; -2,-2; 2,-2; 2,3; 2+a,3; 2+a,-2-a; -2-a, -2-a; -2-a,3];
    case 'peridem_tube'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	s = meshsize;

	contact_radius = 1.74e-04;	% peridem
	space = contact_radius/2;
		% % for rectangular wall, enter internal surface information
		%wall_left = -1 * 1e-3;
		%wall_right = 1 * 1e-3;
		%wall_bottom = -1 * 1e-3;
		%wall_top = 3 * 1e-3;	% fake
		wall_left = -1 * 1e-3-space;
		wall_right = 1 * 1e-3+space;
		wall_bottom = -1 * 1e-3;
		wall_top = 3 * 1e-3;	% fake
	P = [wall_left,wall_top; wall_left,wall_bottom; wall_right,wall_bottom;  wall_right, wall_top; wall_right+s,wall_top; wall_right+s, wall_bottom-s; wall_left-s,wall_bottom-s; wall_left-s, wall_top];
otherwise

end


% Load the geom package

pkg load geometry


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



printf('total nodes: %d \n', length(Pos));



switch geometry
    case {'peridem_wall' , 'peridem_tube'}
	wall_Pos = Pos;
	wall_Vol = Vol;
	wall_T = T;

	printf('Saving matrices: wall_Pos.mat, wall_Vol.mat, wall_T.mat\n')
	% save in matlab-readable format
	save -mat7-binary 'wall_Pos.mat' 'wall_Pos'
	save -mat7-binary 'wall_Vol.mat' 'wall_Vol'
	save -mat7-binary 'wall_T.mat' 'wall_T'


	printf('Saving wall inner dimension: wall_left, right, top, bottom.mat\n')
	% save in matlab-readable format
	save -mat7-binary 'wall_left.mat' 'wall_left'
	save -mat7-binary 'wall_right.mat' 'wall_right'
	save -mat7-binary 'wall_bottom.mat' 'wall_bottom'
	save -mat7-binary 'wall_top.mat' 'wall_top'
	
    otherwise
	printf('Saving matrices: Pos.mat, Vol.mat, T.mat\n')

	% save in matlab-readable format
	save -mat7-binary 'Pos.mat' 'Pos'
	save -mat7-binary 'Vol.mat' 'Vol'
	save -mat7-binary 'T.mat' 'T'
	save -mat7-binary 'delta.mat' 'delta'
	save -mat7-binary 'meshsize.mat' 'meshsize'
	save -mat7-binary 'geometry.mat' 'geometry'

end

disp 'Press key to continue'
pause

% generate the neighbor-list
[NbdArr] = gen_nbdlist(Pos, delta);

[nn, l] = size(NbdArr);

printf('Max neighbors: %d \n', l);
disp 'Press key to continue'
pause
