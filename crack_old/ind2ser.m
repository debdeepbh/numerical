function s=ind2ser(i,j)
global length;
global width;
global dx;
global dy;

nx = length/dx +1;
ny = width/dy+1;

if ( (1 <= i) &&(i <= ny) && (1 <=j) && (j <=nx))
    s = nx*(i-1) + j;
else
    disp('bad value')
end