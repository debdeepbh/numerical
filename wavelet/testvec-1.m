% defining the testvec for fig 23, 24, 31
testvec = zeros(1,512);
for j = 129:256
	n = j-1;
	testvec(j) = sin (((abs(n - 128))^1.7) / 128);
end
for j = 385:448
	n = j-1;
	testvec(j) = sin (((abs(n - 128))^2) / 128);
end



