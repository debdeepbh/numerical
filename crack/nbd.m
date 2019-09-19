function [initial, dx, dy, nx, ny, delta, Xdiff1, Xdiff2, Xdiffnorm] = nbd()
% working with mm
%global Nbd;

length = 0.1; % x axis
width = 0.04;     % y axis


%dx = 0.0005;
%dy = 0.0005;
%
%%delta = 0.002;
%delta = dx*3;
delta =  0.002;
dx = delta/4;
dy = delta/4;


nx = length/dx+1;
ny = width/dy+1;


totalnodes = nx*ny;
%Nbd=zeros(totalNodes,48);   %suggested by matlab for speed. Obviously 2*delta/dx + 2*delta/dy -1 is enough.

tic

for i=1:ny
    for j=1:nx
        iLow = max(1,i-delta/dy);
        iHigh = min(i+delta/dy, ny);
        jLow = max(1, j-delta/dx);
        jHigh = min(j+delta/dx,nx);
        count = 1;
        for k = iLow : iHigh
            for l = jLow : jHigh
               %if ( (((j-l)*dx)^2 + ((k-i)*dy)^2) <= (delta^2) ) %if scaled (k,l) is in delta nbd of scaled (i,j)
               if ( sqrt(((j-l)*dx)^2 + ((k-i)*dy)^2) <= (delta) ) %if scaled (k,l) is in delta nbd of scaled (i,j)
                   if k ~=i || l ~=j
                       if (intersects([0 width/2],[length/2 width/2],ind2pos(i,j,dx,dy),ind2pos(k,l,dx, dy))== 0)
                            Nbd(ind2ser(i,j,nx, ny),count)= (k-1)*nx + l;
                            count=count+1;                           
                       end
                   end
               end
            end
        end
    end
end

Xdiff1 = zeros(size(Nbd))+1;    %replacing unused values with 1
Xdiff2 = zeros(size(Nbd))+1;
Xdiffnorm = zeros(size(Nbd)) +1 ;   %Avoid NaN (in projFA2.m, dividing by this value)

for i=1:totalnodes
    for j = find(Nbd(i,:))
        nbdser = Nbd(i,j);
        x = ser2pos(i,dx,dy, nx, ny);
        xcap = ser2pos(nbdser, dx, dy, nx, ny);
        xdiff = xcap - x;
        xdiffnorm = norm(xdiff);
        Xdiff1(i,j) = xdiff(1);
        Xdiff2(i,j) = xdiff(2);
        Xdiffnorm(i,j) = xdiffnorm;
    end
end
initial = Nbd;

toc
        
        
    

nx
ny
