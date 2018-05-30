% defining the testvec
testvec = zeros(1,512);
for j = 32:95
	testvec(j+1) =	1;
end
for j = 132:259
	testvec(j+1) =	2;
end
for j = 325:415
	testvec(j+1) =	4;
end
