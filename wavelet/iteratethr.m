% iterate over values to get the optimum
ilim = 1:50:1000;
thri = zeros(1,length(ilim));
for i=ilim
	th= testfwd(3, 'd6', 4, [i 10 10 10 500],2);
	thri(i) = th(1);
end
plot(ilim,thri)
