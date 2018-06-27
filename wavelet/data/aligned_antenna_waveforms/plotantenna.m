an = load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g125275173_taxis.txt')(:,2);
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g225275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g325275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g425275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g525275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g625275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g725275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g825275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g925275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g1025275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g1125275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g1225275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g1325275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g1425275173_aligned.txt')];
an = [an load('/home/debdeep/numerical/wavelet/data/aligned_antenna_waveforms/g1525275173_aligned.txt')];

an = an';

%delay = 

size(an)

for i=1:15
	subplot(5,3,i);
	plot(an(i,:));
	xlim([1 260]);
end
