function ratio = damage_test(initial, final, nx, ny)
ratio = sum((final >0),2) ./ sum((initial >0),2); %sum(A,2) gives row sum

damage = 1 - ratio;

A = zeros(ny,nx);
for i=1:ny
    for j=1:nx
	s=ind2ser(i,j,nx,ny);
        A(i,j) = damage(s);
    end
end

imagesc(A)
