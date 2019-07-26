function savenewpos(Nbd, uvec, i, dx, dy, nx, ny, counter, f, filestr)

dotsize = 10;

space = 0.001;
xlimit = [0.1-space 0.2+space];
ylimit = [0-space 0.04+space];

ss = size(uvec);

longVec = (1:ss(1))';

scaling = 10;

XCoord = ser2pos1(longVec, dx, nx) + scaling * uvec(:,1);
YCoord = ser2pos2(longVec, dy, nx) + scaling * uvec(:,2);

%% turn off visible figure
%set(0, 'DefaultFigureVisible', 'off');
%f = figure('visible', 'off');
%set(0, 'currentfigure', f);

filename = strcat('img/', filestr, sprintf('pos%d_%03d.png', i, counter));
if i==3
	Unorm = sqrt(uvec(:,1).^2 + uvec(:,2).^2);
	scatter(XCoord, YCoord, dotsize, Unorm, 'filled')
else
	scatter(XCoord, YCoord, dotsize, uvec(:,i), 'filled')
end
colormap summer;
colorbar;
axis([xlimit ylimit])

saveas(f, filename);
