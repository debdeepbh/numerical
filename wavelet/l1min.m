%% testing the solution
clear all
testvec_7
%for ssize=60:
global samplesize = 200 ;
% .' computes nonconjugate transpopse as opposed to '
N = length(testvec);
%pick a random permulation
global subset = randperm(N, samplesize);
integ = 0:N-1;
global v = 1 - exp(-2*pi*i*integ/N);
%rgiven = (v.*fft(testvec)).';
rgiven = (fft(testvec)).';
global fgiven = rgiven(subset);
%global fgiven = (fft(testvec)).';
global fMat =  dftM(N)(subset,:);

function r = g (x)
	global fgiven;
	global fMat;
	global subset;
	%plot(given)
	%% .' computes nonconjugate transpopse as opposed to '
	r = fMat*x - fgiven;
%  r = [ sumsq(x)-10;
%        x(2)*x(3)-5*x(4)*x(5);
%        x(1)^3+x(2)^3+1 ];
endfunction

function obj = phi (x)
	% total variation norm
	obj = norm(x - circshift(x, [1,0]), 1);
	%% l^2 and TV norm
	%obj = norm(x,2) + norm(x - circshift(x, [1,0]), 1);
	%% l^1 norm
	%obj = norm(x,1);
	%obj = norm(x,2);
%  obj = exp (prod (x)) - 0.5*(x(1)^3+x(2)^3+1)^2;
endfunction

%x0 = [-1.8; 1.7; 1.9; -0.8; -0.8];
%x0 = zeros(size(testvec))'+1;
%x0 = randn(size(testvec))';

x0 = zeros(size(testvec))';
%% not correct since subset is on the fourier side
x0(subset) = (testvec' + (randn(size(testvec))*0.1)')(subset);
%x0 = testvec';

tic
%[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, []);
	%[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g);
	[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g);
toc
