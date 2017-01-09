function A = ser2ind(y)
global length;
global width;
global dx;
global dy;

nx = length/dx +1;
ny = width/dy+1;


if ((y <= (nx*ny)) && (1<=y))
    j = mod(y,nx);
    if j ==0
        j = nx;
    end
    
i = floor((y-1)/nx)+1;
A = [i j];
else
     disp('Bad value');
end
