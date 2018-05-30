% plot the coefficients of wavelet transform w, the p-th stage wavelet transform
function plotcoeffs(w,p)
N = length(w);
type = '.';
for i=1:p
	subplot(p+1,1,i);
	gap = 2^i;
	vals = coeff(w,p,i);
	plot(1:gap:N, vals,type)
	axis([1 N -max(abs(vals)) max(abs(vals))]);
end

%plot the last coarsest part, with the same gap as last
subplot(p+1,1,p+1)
vals = coeff(w,p,p+1);
plot(1:gap:N, vals, type)
axis([1 N -max(abs(vals)) max(abs(vals))])

