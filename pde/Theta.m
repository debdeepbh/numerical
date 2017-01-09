% Theta scheme
function theta()
clc
close all
clear all
b = 0.05;
n = 50;
th = 0.5;
nx = n;
nt = n;
ax = -5;
bx = 5;
at = 0;
bt = 5;
x = linspace (ax, bx, nx);
t = linspace (at, bt, nt);

dt = (bt - at)/(nt-1);
h = (bx -ax)/ (nx-1);
m = dt / (h^2)



% Creating a tridiagonal matrix with diagonal b, subdiagonal b, superdiagonal c with dimension n
%  function T = tridiag(a, b, c, n)
% T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);

A = (1+2*th*b*m)*diag(ones(n,1)) + (-1)*th*b*m*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));

B = (1- 2* 2* (1-th)*m*b)*diag(ones(n,1)) + ((1-th)*m*b)*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));


unot = zeros (1, nx);
k = find (x<=-1);
unot(k) = 0;
k = find (x > -1 & x<=1);
unot (k) = 2*cos(pi * x(k) /2);
k = find (x>1);
unot(k) = 0;

axis([-2 2 0 2])
hold all
plot (x, unot)
pause(1)


plot (x, unot)
pause(1)
%hold all


% Boundary Conditions:

v = unot;

% Loop

for i=1:nx 
	b = B*v';
	w = A \ b;
	v = w';
	plot (x, v);
	
	pause(0.01)
end
	




