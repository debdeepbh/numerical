%p pad by zero
function w = padz(z)
	w = zeros(1,2*length(z));
	% starting index of the actual vector
	% place it in the middle
	start = floor((length(w) - length(z))/2);
	w(start:start+length(z)-1) = z;

