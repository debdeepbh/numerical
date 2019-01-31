% generate the DFT matrix of order n
function M = dftM(N)
W = exp(-2*pi*i/N);
for j = 1:N
	for k =1:N
		M(j,k) = W^((j-1)*(k-1));
	end
end

