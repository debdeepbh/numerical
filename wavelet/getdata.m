% get the wais data from the text files
% use it like wex(3,:), K, noised, noiseund etc
% each column is one data
% after dropping the time data


% the impulse response to be used
K = load('data/notches_260_0_0.txt');
K = K(:,2)';

% get the summed data

i = 1;
wex(i,:,:) = load('data/summed/HpolD21203158.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203755.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203785.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203830.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203866.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203883.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203930.txt')(1:1560,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203958.txt')(1:1560,:); i = i+1;

% drop the timestamps, work unitless
wex = wex(:,:,2);
% now each row is one data

% get the noise data
% un-deconvolved noise
noised = load('data/noise/Deconvolved_HpolD36661151.txt');
noised = noised(:,2)';

noiseund = load('data/noise/observed_HpolC36661151.txt');
noiseund = noiseund(:,2)';
noise = noiseund(320:320+1023);

% get the old theoretical data
testvec_5;	% the spectrum of this one overlaps with that of K
testvecz = testvec;
testconvz = conv(testvecz, K)(1:1024);
testnoisez = randn([1 length(testconvz)])*5;
testz = testconvz + testnoisez;


% get the theoretical data
testvec_6;	% the spectrum of this one overlaps with that of K
testconv = conv(testvec, K)(1:1024);
testnoise = randn([1 length(testconv)])*5;
testy = testconv + testnoise;


% get the wais data manuel gave me in the end

wp1 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203755.txt');
wp2 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203785.txt');
wp3 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203930.txt');
wp4 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203958.txt');
wp5 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29991096.txt');
wp6 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29991181.txt');
wp7 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29991735.txt');
wp8 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29992146.txt');

wpx = [wp1(1:1560,2) wp2(1:1560,2) wp3(1:1560,2) wp4(1:1560,2) wp5(1:1560,2) wp6(1:1560,2) wp7(1:1560,2) wp8(1:1560,2)];

wpx = wpx';
  
% for beamer presentation
z = wpx(4,:)(300:300+1023);
