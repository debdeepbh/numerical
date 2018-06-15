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



% get the theoretical data
testvec_5;
testconv = conv(testvec, K);
testy = testconv + randn([1 length(testconv)])*5;
