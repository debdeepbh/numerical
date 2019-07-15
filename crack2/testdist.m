% test if the neighbors are really within delta distance from the origin
function out = testdist(Nbd)
ss = size(Nbd);

	dist1 = ser2pos1(Nbd) - ser2pos1( (1:ss(1))');
	dist2 = ser2pos2(Nbd) - ser2pos2( (1:ss(1))');

	dist = sqrt(dist1.^2 + dist2.^2);
	
	out = max(max(dist.*(~~Nbd)))

	maxnorm = 0;
	for i=1:ss(1)
		for j = 1:ss(2)
			if Nbd(i,j) > 0
				maxnorm = max(norm(ser2pos(i) - ser2pos(Nbd(i,j))), maxnorm);
			end
		end
	end
