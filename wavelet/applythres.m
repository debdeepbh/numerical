% apply threshold on the p-th stage transform
% thr should be a vector of length p+1
function w = applythres(wo, rule, p, thr)

% don't want to replace the original
w = wo;
N = length(w);

% check
if length(thr) ~= p+1
	disp('threshold vector length incorrect');
end


% staring index
starting = 1;
len = N/2;
for k=1:p
	ending = starting+len-1;
	for l=starting:ending

		switch rule
		case 'hard'
			if( w(l) < thr(k))
				w(l) = 0;
			end
		case 'soft'
			if( w(l) > thr(k))
				w(l) = sign(w(l)) * (abs(w(l)) - thr(k));
			end
		end

	end
	% for the next loop
	len = len/2;
	starting = ending + 1;
end

% for the coarsest level Phi_{-p,k}
for l=starting:N
	if( w(l) < thr(p+1))
		w(l) = 0;
	end
end






