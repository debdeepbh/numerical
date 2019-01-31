
% defining the testvec
testvec = zeros(1,512);
for j = 32:32+4
	testvec(j+1) =	1;
end
for j = 132:132+9
	testvec(j+1) =	2;
end
for j = 325:325+14
	testvec(j+1) =	4;
end
