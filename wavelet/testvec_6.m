% defining the testvec for fig 33 
testvec = zeros(1,512);
for j = 0:511 
       testvec(j+1) = (j - 256)*exp(-(j-256)^2/32);
end



