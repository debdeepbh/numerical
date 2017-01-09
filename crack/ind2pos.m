function C = ind2pos(i,j)
global dx;
global dy;
x = (j-1)*dx;
y = (i-1)*dy;
C = [x y];
