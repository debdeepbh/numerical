function nbdmap(Nbd)
global length;
global width;
global dx;
global dy;

nx = length/dx +1;
ny = width/dy+1;

A = zeros(ny,nx);
for i=1:ny
    for j=1:nx
        serial = ind2ser(i,j);
        row = Nbd(serial,:);
        value = nnz(row);    % number of neighbours
        A(i,j) = value;
    end
end
imagesc(A);

