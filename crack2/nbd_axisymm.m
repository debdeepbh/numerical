function Nbd = nbd_axisymm()
% working with mm
global length;
global width;
global dx;
global dy;
global delta;

global nx;
global ny;

global Xdiffnorm;
global Xdiff1;
global Xdiff2;

global totalnodes;

%global Nbd;

length = 0.1; % x axis
width = 0.04;     % y axis
dx = 0.0005;
dy = 0.0005;
delta = 4 * dx; 	% 0.002;


nx = length/dx+1;
ny = width/dy+1;
totalnodes = nx*ny;
%Nbd=zeros(totalNodes,48);   %suggested by matlab for speed. Obviously 2*delta/dx + 2*delta/dy -1 is enough.

tic

for i=1:ny	% i is the index associated with the y coordianate
    for j=1:nx	% j is for x coordinate
        iLow = max(1, i-delta/dy);
        iHigh = min(i+delta/dy, ny);
        jLow = max(1, j-delta/dx);
        jHigh = min(j+delta/dx, nx);
        count = 1;
        for k = iLow : iHigh
            for l = jLow : jHigh
               if (sqrt(((j-l)*dx)^2 + ((k-i)*dy)^2) < delta) %if scaled (k,l) is in delta nbd of scaled (i,j)

                   if k ~=i || l ~=j
                      % if (intersects([0 width/2],[length/2 width/2],ind2pos(i,j),ind2pos(k,l))== 0)
                            Nbd(ind2ser(i,j),count)= ind2ser(k,l); 
                            count=count+1;                           
                       %end
                   end
               end
            end
        end
    end
end

%Xdiff1 = zeros(size(Nbd))+1;    %replacing unused values with 1
%Xdiff2 = zeros(size(Nbd))+1;
%Xdiffnorm = zeros(size(Nbd)) +1 ;   %Avoid NaN (in projFA2.m, dividing by this value)
%
%for i=1:totalnodes
%    for j = find(Nbd(i,:))
%        nbdser = Nbd(i,j);
%        x = ser2pos(i);
%        xcap = ser2pos(nbdser);
%        xdiff = xcap - x;
%        xdiffnorm = norm(xdiff);
%        Xdiff1(i,j) = xdiff(1);
%        Xdiff2(i,j) = xdiff(2);
%        Xdiffnorm(i,j) = xdiffnorm;
%    end
%end


toc
        
        
    

