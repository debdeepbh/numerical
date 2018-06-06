%plot the unsummed, unconvolved signals seen by a particular antenna
a(1,:,:) = load('/home/debdeep/gdrive/anita/wais_pulse/unsummed/antennaGraphs_21203883/ant1.csv');
a(2,:,:) = load('/home/debdeep/gdrive/anita/wais_pulse/unsummed/antennaGraphs_21203883/ant2.csv');
a(3,:,:) = load('/home/debdeep/gdrive/anita/wais_pulse/unsummed/antennaGraphs_21203883/ant3.csv');
a(4,:,:) = load('/home/debdeep/gdrive/anita/wais_pulse/unsummed/antennaGraphs_21203883/ant4.csv');

for i = 1:4
    subplot(4,1,i);
    plot(a(i,:,1),a(i,:,2));
    axis([5 100 -60 60]);
end


