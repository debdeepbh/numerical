% defining the testvec for fig 33
testvec = zeros(1,512);
for j = 0:511
	testvec(j+1) = sin(j^(1.5)/64);
end



