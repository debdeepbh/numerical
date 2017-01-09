% FTFS
function Forward()
clc
close all
clear all
c = -1;
nx = 101;
nt = 30;
ax = -2;
bx = 2;
at = 0;
bt = 1;
x = linspace (ax, bx, nx);
t = linspace (at, bt, nt);
dt = (bt-at)/(nt-1);

m = c*((bt-at)*(nx-1))/((bx-ax)*(nt-1));


unot = zeros (1, nx);
k = find (x<=0);
unot(k) = 1;
k = find (x > 0 & x<=1);
unot(k) = 1 - x(k);
k = find (x>1);
unot(k) = 0;
plot (x, unot)
pause(1)
%hold all

v = zeros (nt , nx);
v(1,:) = unot;

% Boundary Conditions:

%v(:, 1) = 1;

for n = 1:nt-1
	for i = 1:nx-1
	v(n+1,i)= ((1+m)*v(n,i))-m*v(n,i+1);
	end 

end

for j = 2:nt
plot (x, v(j,:))
pause(0.1)
%hold all
end


