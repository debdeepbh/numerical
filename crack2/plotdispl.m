function plotdispl(uvec, i)

dotsize = 20;

ss = size(uvec);

longVec = (1:ss(1))';

XCoord = ser2pos1(longVec);
YCoord = ser2pos2(longVec);


if i==3
	Unorm = sqrt(uvec(:,1).^2 + uvec(:,2).^2);
	scatter(XCoord, YCoord, dotsize, Unorm, 'filled')
else
	scatter(XCoord, YCoord, dotsize, uvec(:,i), 'filled')
end

colormap summer;
	colorbar
