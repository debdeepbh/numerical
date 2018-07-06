% prepare a signal with necessary data
function [ax, aximp] = prepsig(num)


datapath = '/home/debdeep/numerical/wavelet/data/waveforms/info_';

pathname = '/home/debdeep/numerical/wavelet/data/notches_Debdeep/';

% naming of the signals: ax(i,:)
% naeing of the impulse response: aximp(i,:)
ax = load(strcat(datapath,num2str(num),'_sig'))';
antlist = load(strcat(datapath,num2str(num),'_ant'));

aximp = zeros(length(antlist),1024);
for i = 1:length(antlist)
	aximp(i,:) = load(strcat(pathname,getantenna(antlist(i))))(:,2)';
end
