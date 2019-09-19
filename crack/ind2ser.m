function s=ind2ser(i,j,nx,ny)

if ( (1 <= i) &&(i <= ny) && (1 <=j) && (j <=nx))
   s = nx*(i-1) + j;
else
	nx, ny
	i,j
    disp('bad value')
    return
end
