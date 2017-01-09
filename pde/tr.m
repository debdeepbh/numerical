function tr()
clc 
close all 
clear all
% solving the problem
% u_t + u_x = 0


% Mesh for spacs
a = -2;
b = 2;
nnode = 100;
x = linspace(a, b, nnode);
% x = a:0.5:b
h = (b-a)/(nnode-1);

% time parameter
Tf = 1;
CFL = 0.9;
c = -1;
dt = CFL*h/abs(c);

% Initial condition

u0 = zeros(nnode,1);

for i = 1:nnode
	if (x(i) <= 1)
		u0(i) = 1;
	elseif(x(i) <= 1)
		u0(i) = 1 - x(i);
	end
end

%Solve with forward time forward space FTFS Algorithm


un = u0;
unp1 = u0;

plot(x,u0)
	title('solutions')
	xlabel('x')
	ylabel('Solution')
	pause(2)

t = 0;
lambda = c*dt/h;
while( t <= Tf)
	unp1(1:nnode-1) = un(1:nnode-1) - lambda*(un(2:nnode) - un(1:nnode-1)); 

	plot(x,unp1)
	title('solutions')
	xlabel('x')
	ylabel('Solution')
	pause(0.1)

	t = t + dt;

	un = unp1;
end












