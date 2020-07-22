
% generate the mesh

%script

clear all
close all



geometry='E'
%geometry='heart'
%geometry='leaf'
%geometry='peridem'
%geometry='pacman'
%geometry='peridem_wall'
%geometry='peridem_tube'
%geometry='peridem_tube_onesided'
%geometry = 'peridem_w_prenotch'
%geometry = 'peridem_triangle'
%geometry = 'peridem_floor'
%geometry='unitcircle'
%geometry='sodalime'
%geometry='big_box'
%geometry='pile_box'
%geometry='pile_box_narrow'
%geometry='pile_ball'
%geometry='moving_wall'


draw_text = 1;

switch geometry
case	'E'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
P = 1/2* [0.02,1.4; 0.02,-1.4; 2.02,-1.4; 2.02,-1.0; 0.42,-1.0; 0.42,-0.2; 1.22,-0.2; 1.22,0.2; 0.42,0.2; 0.42,1.0; 2.02,1.0; 2.02,1.4] * 1e-3;
case 'heart'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
    P = 1/3 * (1e-3 .* [3.2215934,1.0041422; 2.8215933,1.8041421; 2.0215933,2.204142; 1.2215934,2.204142; 0.42159337,1.8041421; 0.021593383,1.0041422; 0.42159337,-0.19585787; 1.6215934,-0.99585783; 2.8215933,-1.7958579; 3.2215934,-2.1958578; 3.6215935,-1.7958579; 6.0215936,-0.19585787; 6.421593,1.0041422; 6.0215936,1.8041421; 5.2215934,2.204142; 4.421593,2.204142; 3.6215935,1.8041421] - [0.0032, 0]);

case	'leaf'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
    P = 1e-3 *[2.8624094,3.1853633; 2.0624094,1.1853632; 1.2624093,1.5853633; 1.6624093,-0.014636788; 0.4624093,0.38536322; 1.6624093,-1.2146368; 0.0624093,-1.2146368; 2.8624094,-3.2146368; 5.6624093,-1.2146368; 4.0624094,-1.2146368; 4.862409,0.38536322; 3.6624093,-0.014636788; 4.4624095,1.5853633; 3.6624093,1.1853632]
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
    case 'big_box'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	s = meshsize;
		% for rectangular wall
		wall_left = -5 * 1e-3;
		wall_right = 5 * 1e-3;
		wall_top = 5 * 1e-3;	% fake
		wall_bottom = -5 * 1e-3;
	P = [wall_left,wall_top; wall_left,wall_bottom; wall_right,wall_bottom;  wall_right, wall_top; wall_right+s,wall_top; wall_right+s, wall_bottom-s; wall_left-s,wall_bottom-s; wall_left-s, wall_top];

    case 'pile_box_narrow'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	s = meshsize;
		% for rectangular wall
		wall_left = -2 * 1e-3;
		wall_right = 2 * 1e-3;
		wall_top = 5 * 1e-3;	% fake
		wall_bottom = -5 * 1e-3;
	P = [wall_left,wall_top; wall_left,wall_bottom; wall_right,wall_bottom;  wall_right, wall_top; wall_right+s,wall_top; wall_right+s, wall_bottom-s; wall_left-s,wall_bottom-s; wall_left-s, wall_top];

    case 'moving_wall'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
		% for rectangular wall
		moving_wall_left = -2 * 1e-3;
		moving_wall_right = 2 * 1e-3;
		moving_wall_top = 3 * 1e-3;	% fake
		moving_wall_bottom = 2.8 * 1e-3;
	P = [moving_wall_left, moving_wall_bottom; moving_wall_right, moving_wall_bottom; moving_wall_right, moving_wall_top; moving_wall_left, moving_wall_top];

    case 'pile_box'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/(2);
	% thickness
	s = meshsize;
		% for rectangular wall
		wall_left = -4 * 1e-3;
		wall_right = 4 * 1e-3;
		wall_top = 5 * 1e-3;	% fake
		wall_bottom = -5 * 1e-3;
	P = [wall_left,wall_top; wall_left,wall_bottom; wall_right,wall_bottom;  wall_right, wall_top; wall_right+s,wall_top; wall_right+s, wall_bottom-s; wall_left-s,wall_bottom-s; wall_left-s, wall_top];

    case 'pile_ball'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/(2);

	% circle, Caution: the last node should _NOT_ overlap with the first!
	steps = 20;
	angles = linspace(0, 2*pi- 2*pi/steps, steps)';
	P = 1e-3 * [ cos(angles), sin(angles)]; % 1mm size

    case 'peridem_tube'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	s = meshsize;

	contact_radius = 1.74e-04;	% peridem
	space = contact_radius/2;
	%space = contact_radius;
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
    case 'peridem_tube_onesided'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	s = meshsize;

	contact_radius = 1.74e-04;	% peridem
	space = contact_radius/2;
	%space = contact_radius;
		% % for rectangular wall, enter internal surface information
		%wall_left = -1 * 1e-3;
		%wall_right = 1 * 1e-3;
		%wall_bottom = -1 * 1e-3;
		%wall_top = 3 * 1e-3;	% fake
		wall_left = -1 * 1e-3-space;
		wall_right = 2 * 1e-3+space;
		wall_bottom = -1 * 1e-3;
		wall_top = 3 * 1e-3;	% fake
	P = [wall_left,wall_top; wall_left,wall_bottom; wall_right,wall_bottom;  wall_right, wall_top; wall_right+s,wall_top; wall_right+s, wall_bottom-s; wall_left-s,wall_bottom-s; wall_left-s, wall_top];

    case 'peridem_floor'
	delta = 1e-3;	% peridem, all the nodes are neighbors
	meshsize = delta/5;
	% thickness
	s = meshsize;

	contact_radius = 1.74e-04;	% peridem
	space = contact_radius/2;
	%space = contact_radius;
		% % for rectangular wall, enter internal surface information
		%wall_left = -1 * 1e-3;
		%wall_right = 1 * 1e-3;
		%wall_bottom = -1 * 1e-3;
		%wall_top = 3 * 1e-3;	% fake
		wall_left = -2e-3;
		wall_right = 2e-3;
		wall_bottom = -0.25e-3;
		wall_top = 0e-3;
	P = [wall_left, wall_bottom; wall_right, wall_bottom; wall_right, wall_top; wall_left, wall_top];
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
    case {'peridem_wall' , 'peridem_tube', 'peridem_tube_onesided', 'peridem_floor', 'pile_box', 'pile_box_narrow'}
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
    case {'moving_wall'}
	moving_wall_Pos = Pos;
	moving_wall_Vol = Vol;
	moving_wall_T = T;

	printf('Saving matrices: moving_wall_Pos.mat, moving_wall_Vol.mat, moving_wall_T.mat\n')
	% save in matlab-readable format
	save -mat7-binary 'moving_wall_Pos.mat' 'moving_wall_Pos'
	save -mat7-binary 'moving_wall_Vol.mat' 'moving_wall_Vol'
	save -mat7-binary 'moving_wall_T.mat' 'moving_wall_T'


	printf('Saving moving_wall inner dimension: moving_wall_left, right, top, bottom.mat\n')
	% save in matlab-readable format
	save -mat7-binary 'moving_wall_left.mat' 'moving_wall_left'
	save -mat7-binary 'moving_wall_right.mat' 'moving_wall_right'
	save -mat7-binary 'moving_wall_bottom.mat' 'moving_wall_bottom'
	save -mat7-binary 'moving_wall_top.mat' 'moving_wall_top'
	
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
