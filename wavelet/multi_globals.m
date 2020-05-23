%% set up the global variables
	 N = 1024;
% using global variables is a bad idea it seems

	 M = 15;
	 type='d10';
	% type='meyer';
	 p=3;
	 rho=1;
	 method='hard';
	 B = getbasismat(type, p, N);
	 fB = fft(B')';

% running the same script again with one of the global variables changed does not change the value of the variable (if called without global)
	%global M = 15;
	%global type='d10';
	%%global type='meyer';
	%global p=3;
	%global rho=1;
	%global method='hard';
	%global B = getbasismat(type, p, N);
	%global fB = fft(B')';

	


