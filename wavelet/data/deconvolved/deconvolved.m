wp1 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203755.txt');
wp2 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203785.txt');
wp3 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203830.txt');
wp4 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203866.txt');
wp5 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203883.txt');
wp6 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203930.txt');
wp7 = load('/home/debdeep/gdrive/anita/wais_pulse/deconvolved/HpolD21203958.txt');
wp = [wp1(1:1560,1) wp1(1:1560,2) wp2(1:1560,2) wp3(1:1560,2) wp4(1:1560,2) wp5(1:1560,2) wp6(1:1560,2) wp7(1:1560,2)];



figure;
for i=2:8
	subplot(7,1,i);
	plot(wp(:,i));
	grid on;
end