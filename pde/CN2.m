% Crank-Nicolson using Gauss Elimination
function CN2()
clc
close all
clear all
c = 1;
n = 50;

nx = n;
nt = n;
ax = -2;
bx = 2;
at = 0;
bt = 1;
x = linspace (ax, bx, nx);
t = linspace (at, bt, nt);
dt = (bt-at)/(nt-1);

m = c*((bt-at)*(nx-1))/((bx-ax)*(nt-1));



% Creating a tridiagonal matrix with diagonal b, subdiagonal b, superdiagonal c with dimension n
%  function T = tridiag(a, b, c, n)
% T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);

A = diag(ones(n,1)) + (m/4)*diag(ones(n-1,1),1) + (-1)*(m/4)*diag(ones(n-1,1),-1);

unot = zeros (1, nx);
k = find (x<=-1);
unot(k) = 0;
k = find (x > -1 & x<=1);
unot (k) = cos(pi * x(k) /2);
k = find (x>1);
unot(k) = 0;
plot (x, unot)
pause(1)
%hold all


% Boundary Conditions:

v = unot;

% Loop

for i=1:nx 
	b = A'*v';
	w = A \ b;
	v = w';
	plot (x, v);
	
	pause(0.1)
end
	




