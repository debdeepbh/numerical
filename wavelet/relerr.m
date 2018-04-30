% returns the relative error matrix in pth stage wavelet compression, with respect to K, in q norm
function Rel = relerr(z, typelist, maxK, p, q)
% e.g. typelist = {'shan'; 'db5'}
% number of methods to compare
r = length(typelist);
% zero vector 
Rel = zeros(r,maxK);
for K=0:maxK-1
	% here typelist is a cell array. Use stringlist{1} etc
	for j = 1:r
		Rel(j,K+1) = norm(z - compress(z, typelist{j}, p, K), q)/norm(z, q);
	end
end

%%print(tf,'-dpng');
