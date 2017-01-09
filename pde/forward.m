% Solving transport equation using forward space forward time method
function forward()
clear all
close all
nx = 5;
nt = 10;
xa = -5;
xb = 5;
ta = 0;
tb = 10;
dx = (xb - xa)/nx;
dt = (tb - ta)/nt;
lambda = dt/dx
x = linspace(xa, xb, nx);
t = linspace(ta, tb, nt);

unot = zeros(1, nx);
k = find(x<0);
unot(k)= 0;
k = find(0<x);
unot(k) =   x(k);

%plot (x, unot)

hold all
v = zeros(nt, nx);
v(1,:)=unot;
%v(:,nx)= 0;
plot(x,v(1,:))
for n =2:nt
	
	for i =1:nx-1
		v(n,i)=v(n-1,i) - lambda * (v(n-1,i+1) - v(n-1,i));
	end
	v(n,nx)=unot(nx - dt *(n-1));
	%v
	if(n<=10)
	plot(x,v(n,:))
	end
end
%v(2,:);


plot(x,v(nt,:))
axis([-5 5 -10 10])
