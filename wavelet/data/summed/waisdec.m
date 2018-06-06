% deconvolve using the notches_260_0_0 files
waispulse

K = load('notches_260_0_0.txt');
K = K(:,2);
figure;
% have the matrix wp
for i=2:8
	sig = wp(:,i);
	subplot(7,1,i);
	plot(decmore(sig',K','allp'));
end


