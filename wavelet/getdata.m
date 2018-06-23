% get the wais data from the text files
% use it like wex(3,:), K, noised, noiseund etc
% each column is one data
% after dropping the time data


% the impulse response to be used
K = load('data/notches_260_0_0.txt');
K = K(:,2)';

% get the summed data

i = 1;
starting = 1;
ending = 1560;
wex(i,:,:) = load('data/summed/HpolD21203158.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203755.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203785.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203830.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203866.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203883.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203930.txt')(starting:ending,:); i = i+1;
wex(i,:,:) = load('data/summed/HpolD21203958.txt')(starting:ending,:); i = i+1;

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
testzori = [testvecz zeros(1,512)];


% get the theoretical data
testvec_6;	% the spectrum of this one overlaps with that of K
testconv = conv(testvec, K)(1:1024);
testnoise = randn([1 length(testconv)])*5;
testy = testconv + testnoise;
testyori = [testvec zeros(1,512)];


% get the wais data manuel gave me in the end
starting = 320;
ending = 320+1023;
wp1 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203755.txt');
wp2 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203785.txt');
wp3 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203930.txt');
wp4 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC21203958.txt');
wp5 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29991096.txt');
wp6 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29991181.txt');
wp7 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29991735.txt');
wp8 = load('/home/debdeep/gdrive/anita/wais_pulse/debd_stuff/HpolC29992146.txt');

wpx = [wp1(starting:ending,2) wp2(starting:ending,2) wp3(starting:ending,2) wp4(starting:ending,2) wp5(starting:ending,2) wp6(starting:ending,2) wp7(starting:ending,2) wp8(starting:ending,2)];

wpx = wpx';
  
% for beamer presentation
z = wpx(4,:);
