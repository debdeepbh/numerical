%=================================================================================
% Solving the problem
% u_t + c u_x = 0
% u(x,0) = u0(x)
%=================================================================================

% Generating the mesh/grid on the domain [a,b] = [-2,2]
%function dip()

clc
clear all
close all

c = -1;


a = -2;
b = 2;
nnode = 100;
x = linspace(a,b,nnode);
h = (b-a)/(nnode-1);


% Time Parameter
Tf = 3;

% CFL number: |c|*dt/h where dt/h is called the Mesh Speed
% This has to be less than one
CFL = 1;
dt = CFL*h/abs(c);

% Initial condition
u0 = zeros(nnode,1);
for i =1:nnode
	if((x(i) >= -1) & (x(i) <= 0))
		u0(i) = 1 + x(i);
	elseif ((x(i) <= 1)&(x(i) >=0))
		u0(i) = 1 - x(i);
	end
end

% Solve with forward time forward space algorithm

un = u0;	% current time
unpl = u0;	% future time

% plot once to see
	plot(x,u0)
	title('initial value')
	xlabel('x')
	ylabel('solution')
	pause(2);


t = 0;
lambda = c*dt/h;
while (t  <= Tf)
	unpl(1:nnode-1) = un(1:nnode-1) - lambda * (un(2:nnode) - un(1:nnode -1));
	plot(x,unpl)
	title('solutions')
	xlabel('x')
	ylabel('solution')
	pause (0.1);
	t = t + dt;
	
	un = unpl;

end
	
	
