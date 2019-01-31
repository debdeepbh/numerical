%% given length, total levels, produces an array of indices corresponding to the qth level
function [startInd, lastInd] = getindex(N, p, q)
if q == p+1
	lastInd = N;
	startInd = N - N/2^p +1;
else

	lastInd =0;
	prevLen = N;
	for iter=1:q
		startInd = lastInd + 1;
		levelLength = prevLen/2;
		lastInd = startInd + levelLength -1;
		prevLen = levelLength;
	end
end





