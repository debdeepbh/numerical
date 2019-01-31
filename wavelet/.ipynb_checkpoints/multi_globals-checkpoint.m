%% set up the global variables
	global N = 1024;
	global M = 15;
	global type='meyer';
	global p=3;
	global rho=1;
	global method='hard';
	global B = getbasismat(type, p, N);
	global fB = fft(B')';

	


