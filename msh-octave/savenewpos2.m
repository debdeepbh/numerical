function savenewpos2(uvec, Pos, i, counter, f, filestr)
% Save current positions to a file. n=#nodes
% Input:
%	uvec: displacement vector, nx2
%	Pos: position vector of the nodes in the reference configuration, nx2
%	i: which quantity to represent in the color, i=1,2 are x,y coords, 3=norm
%	counter: plot counter
%	filestr: string prefix for the saved filenames
%
% Output:
%	:

dotsize = 10;


%scaling = 5;
scaling = 10;

CurrPos = Pos + scaling * uvec;

%% turn off visible figure
%set(0, 'DefaultFigureVisible', 'off');
%f = figure('visible', 'off');
%set(0, 'currentfigure', f);

filename = strcat('img/', filestr, sprintf('pos%d_%03d.png', i, counter));
if i==3
	u_norm = sqrt(uvec(:,1).^2 + uvec(:,2).^2);

	scatter(CurrPos(:,1), CurrPos(:,2), dotsize, u_norm, 'filled')
else
	scatter(CurrPos(:,1), CurrPos(:,2), dotsize, uvec(:,i), 'filled')
end
colormap summer;
colorbar;

space = 0.001;
%xmin = min(Pos(:,1)) - space;
%ymin = min(Pos(:,2)) - space;
%xmax = max(Pos(:,1)) + space;
%ymax = max(Pos(:,2)) + space;
%xlimit = [xmin xmax];
%ylimit = [ymin ymax];
%axis([xlimit ylimit])

axis equal

saveas(f, filename);
