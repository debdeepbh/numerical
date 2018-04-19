% defining the testvec for28, 29, 30
testvec = zeros(1,512);
for j = 0:63
	testvec(j+1) = 1 - j/64;
end
for j = 256:319
	testvec(j+1) = 5 - j/64;
end
for j = 64:255
	testvec(j+1) =0;
end
for j = 320:511
	testvec(j+1) =0;
end



