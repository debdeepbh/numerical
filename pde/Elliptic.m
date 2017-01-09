% Solution to elliptic equation Laplacian u = f, u = g on bdry
function Elliptic()
clc
close all
clear all
M = 10;
N = M*M;

nx = M;
nt = M;
ax = 0;
bx = 1;
at = 0;
bt = 1;
x = linspace (ax, bx, nx);
t = linspace (at, bt, nt);
dt = (bt-at)/(nt-1);
h = (bx -ax)/(nx -1);

f = (-1)*ones(N,N);
g = 0;

A = zeros(N,N);
b = zeros(N,1);

for n=1:N
	j = floor((n-1)/M) + 1;
	i = rem((n-1), M) + 1;
	if(i==1 | i==M | j==1 | j==M)
		A(n,n) = 1;
		b(n,1) = g;
	else
		A(n,n) = 4;
		A(n,n-1) = -1;
		A(n,n+1) = -1;
		A(n, (j-2)*M + i) = -1;
		A(n, j*M + i) = -1;
		
		b(n,1) = (-1)*h*h*f(i,j);
	end
end

w = A \ b;
v = zeros(M,M);
for n = 1:N
	j = floor((n-1)/M) + 1;
	i = rem((n-1), M) + 1;
	v(i,j) = w(n,1);	
end

surf(x,t,v'); 	% Had to take transpose of v since:
		% x = x(i), y = y(j); r = r(i,j)
		% We have to plot the tuples (x(i),y(j),r(i,j)) in R^3
		% plot3(x,y,r') or surf(x,y,r') gives the desired result

