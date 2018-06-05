% get the projection of a vector onto the p-th level wavelet basis
% p can be from 0 to log2(N)
function prj = proj(z, type, p)

% get the wavelet transform
w = wtrans(z, type, p);

N = length(w);

% staring index
starting = 1;
len = N/2;
for k=1:p
	ending = starting+len-1;
	w(starting:ending) = 0;
	% for the next loop
	len = len/2;
	starting = ending + 1;
end

prj = iwtrans(w, type, p);





