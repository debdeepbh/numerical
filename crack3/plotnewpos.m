function plotnewpos(Nbd, uvec, i, dx, dy, nx, ny)

dotsize = 100;

ss = size(uvec);

longVec = (1:ss(1))';

scaling = 10;

XCoord = ser2pos1(longVec, dx, nx) + scaling * uvec(:,1);
YCoord = ser2pos2(longVec, dy, nx) + scaling * uvec(:,2);


if i==3
	Unorm = sqrt(uvec(:,1).^2 + uvec(:,2).^2);
	scatter(XCoord, YCoord, dotsize, Unorm, 'filled')
else
	scatter(XCoord, YCoord, dotsize, uvec(:,i), 'filled')
end

colormap summer;
	colorbar

