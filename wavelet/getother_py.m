% get v from my notation using the values listed in python library for v in their notation
% http://wavelets.pybytes.com/wavelet/bior1.3/
% Where the zeros added in the beginning to construct 
function v = getother_py(u) 
N = length(u);
% v(k) = (-1)^{k-1} u_bar(1-k)

% v(k) = (-1)^(2k-3) P(3-k) where P = [zeros Py_v]
for k=0:N-1
	% mod keep the value in the range of 0 and N-1
	% +1 shifts the vector so that u is defined
	% place the result in v(k+1)
	v(k+1) = (-1)^(2*k-3)*conj(u(mod(3-k,N) +1));
end

