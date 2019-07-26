function [initial, dx, dy, nx, ny, delta] = nbd_test()
% working with mm

%% Our y axis extends downward, x axis goes to the right as usual

length = 0.1; % x axis
width = 0.04;     % y axis


%dx = 0.0005;
%dy = 0.0005;
%delta = dx*4;
%% delta = 0.002;

delta =  0.002;
dx = delta/4;
dy = delta/4;


nx = length/dx+1;
ny = width/dy+1;


totalnodes = nx*ny;
%Nbd=zeros(totalNodes,48);   %suggested by matlab for speed. Obviously 2*delta/dx + 2*delta/dy -1 is enough.

midserial_i =  (ny -1)/2+1;
midserial_j =  (nx -1)/2+1;

midpoint_posx = indj2pos1(midserial_j, dx);
midpoint_posy = indi2pos2(midserial_i, dy);

bottom_left_posx = indj2pos1(1, dx);
bottom_left_posy = indi2pos2(ny, dy);

disp('Bottom left in upside-down coordinate system')
[bottom_left_posx bottom_left_posy]
disp('Midpoints')
[midpoint_posx midpoint_posy]


tic



for i=1:ny
	fprintf(repmat('\b',1,100));
	fprintf('Generating Mesh: %d %% ',floor(i/ny*100));

    for j=1:nx
        iLow = max(1,i-delta/dy);
        iHigh = min(i+delta/dy, ny);
        jLow = max(1, j-delta/dx);
        jHigh = min(j+delta/dx,nx);
        count = 1;
        for k = iLow : iHigh
            for l = jLow : jHigh
               if ( sqrt(((j-l)*dx)^2 + ((k-i)*dy)^2) <= (delta) ) %if scaled (k,l) is in delta nbd of scaled (i,j)
                   if k ~=i || l ~=j
%                       if (intersects([bottom_left_posx midpoint_posy], [midpoint_posx midpoint_posy], [indj2pos1(j,dx) indi2pos2(i,dy)], [indj2pos1(l,dx) indi2pos2(k,dy)])== 0)
                            Nbd(ind2ser(i,j,nx, ny),count)= (k-1)*nx + l;
                            count=count+1;                           
%                       end
                   end
               end
            end
        end
    end
end

initial = Nbd;

toc
        
        
    

nx, ny
