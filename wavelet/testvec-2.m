% defining the testvec
testvec = zeros(1,512);
for j = 32:95
       testvec(j+1) = (j - 256)*exp(-(j-256)^2/512);
end



