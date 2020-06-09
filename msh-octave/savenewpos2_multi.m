function savenewpos2_multi(total_particles, uvec, Pos, j, counter, f, filestr)
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


scaling = 5;
%scaling = 10;

CurrPos = Pos + scaling * uvec;

filename = strcat('img/', filestr, sprintf('pos%d_%03d.png', j, counter));
switch j
    case 3
	for i = 1:total_particles
	    u_norm = sqrt(uvec(:,1, i).^2 + uvec(:,2, i).^2);
	    scatter(CurrPos(:,1, i), CurrPos(:,2, i), dotsize, u_norm, 'filled')

	    hold on
	end
	hold off
    otherwise
	for i = 1:total_particles
	    scatter(CurrPos(:,1,i), CurrPos(:,2,i), dotsize, uvec(:,j, i), 'filled')
	    hold on
	end
	hold off
	
end

%colormap summer;
%colormap winter;
%colormap lines;
%colormap jet;
colormap hsv;
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
