function savenewpos2_multi_with_wall(total_particles, CurrPos, wall_CurrPos, Quantity, counter, f, filestr, contact_radius, time)
% Save current positions to a file. n=#nodes
% Input:
%	CurrPos: position vector, nx2xtotal_particle
%	Quantity: quantity to plot, nx1xtotal_particle
%	counter: plot counter
%	filestr: string prefix for the saved filenames
%
% Output:
%	:


dotsize = 10;

filename = strcat('img/', filestr, sprintf('pos%d_%03d.png', j, counter));

	for i = 1:total_particles
	    scatter(CurrPos(:,1, i), CurrPos(:,2, i), dotsize, Quantity(:,:,i), 'filled')
	    hold on
	end


	% for wall, plot the position in color, for now
	scatter(wall_CurrPos(:,1), wall_CurrPos(:,2), dotsize, wall_CurrPos(:,2), 'filled')

	% plot contact radius
	% around this particle
	particle_index = 1;
	% centered at this node
	top_node = 6;
	bottom_node = 16;
	theta = 0:0.1:2*pi;
	plot(CurrPos(bottom_node,1, particle_index) + contact_radius * cos(theta), CurrPos(bottom_node,2, particle_index) + contact_radius * sin(theta));

	hold off

colormap summer;
%colormap winter;
%colormap lines;
%colormap jet;
%colormap hsv;
colorbar;


title(strcat('time = ', num2str(time), ' ms'))

space = 0.001;
%xmin = min(Pos(:,1)) - space;
%ymin = min(Pos(:,2)) - space;
%xmax = max(Pos(:,1)) + space;
%ymax = max(Pos(:,2)) + space;
%xlimit = [xmin xmax];
%ylimit = [ymin ymax];
%axis([xlimit ylimit])

xlim([-2.5 2.5]* 1e-3)
ylim([-3.5 5]* 1e-3)

axis equal
grid on

saveas(f, filename);
