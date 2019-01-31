% Coherence map
% inputPos should be a Mx2 matrix with (i,:) containting the localtion of the ith antenna
function sumall = beam(input, inputPos)
[M, N] = size(input);

ind1 = 1;
ind2 = 1;

for r1 = -1:0.1:1 
	for r2 = -1:0.1:1
%for r1 = -pi:0.2:pi 
%	for r2 = -pi/4:0.1:pi/4
% r = [cos(r2)*cos(r1), cos(r2)*sin(r1)];
r = [ r1, r2 ];

		sumall(ind1, ind2) = 0;
		for i=1:M
			for j=1:i-1
				sumall(ind1, ind2) = sumall(ind1, ind2) + crossc(input(i,:), input(j,:), inputPos(i,:), inputPos(j,:), r);
			end
		end

		ind2 = ind2 + 1;
	end
	ind1 = ind1 + 1;
end
