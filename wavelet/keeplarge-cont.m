% zero out all but highest k elements in a vector
function z = keeplarge(w,k)
N = length(w);
if(k!=0)
	% sort in *ascending* order wrt absolute values
	h = sort(abs(w));
	% find the indices of all (k-2) largest elements
	q1 = find(abs(w) < h(N+1-k));
	% find indices of only the (k-1)-th largest element
	q2 = find(abs(w) == h(N+1-k));
	% in case there are multiple (K-1)-th largest elements, choose some until the total length
	q3 = q2(1:N-k-length(q1));
	q = [q1, q3];

	z = w;
	% use (linear) continuous cutoff instead
	for j=q
		%z(q) = 0;
		weight = 1 - (h(N+1-k) - abs(z(j)))/h(N+1-k);
		z(j) = z(j)* weight/9;
	end
else
	z = zeros(1,N);
end


