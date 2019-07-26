function ratio = damage(initial, final)
ratio = sum((final >0),2) ./ sum((initial >0),2); %sum(A,2) gives row sum

global length;
global width;
global dx;
global dy;

nx = length/dx +1;
ny = width/dy+1;

A = zeros(ny,nx);
for i=1:ny
    for j=1:nx
        A(i,j) = ratio(ind2ser(i,j));
    end
end
imagesc(A)
