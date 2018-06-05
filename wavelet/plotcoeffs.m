% plot the coefficients of wavelet transform w, the p-th stage wavelet transform
function plotcoeffs(w,p)
N = length(w);
type = '.';

% 
figure;


for i=1:p
	subplot(p+1,1,i);
	gap = 2^i;
	vals = coeff(w,p,i);
	stem(1:gap:N, vals,type)
	axis([1 N -max(abs(vals)) max(abs(vals))]);
	% to remove the x-axis labels, markers
	% set(gca, 'XTickLabel', [],'XTick',[]);
end

%plot the last coarsest part, with the same gap as last
subplot(p+1,1,p+1)
vals = coeff(w,p,p+1);
stem(1:gap:N, vals, type)
axis([1 N -max(abs(vals)) max(abs(vals))])

