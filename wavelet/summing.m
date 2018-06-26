% add noise to each level and then sum it up
sigma = 5;
nlen = 15;

% for wavelet
type = 'meyer';
p = 5;
rho = 1;
% sc_a = 1;	% optimal scaling to be computed once
method = 'soft';

N = length(testconv);

noimat = testymat = decf = decw = deca =  zeros(nlen,N);
testysum = noisum = decfsum = decwsum = decasum = zeros(1,N);

for i=1:nlen
	noimat(i,:) = randn([ 1 1024])*sigma;
	testymat(i,:) = testconv + noimat(i,:);

	%% get optimum scaling for the first time
	%if (i==1)
	%	sc_a = getoptsc(testymat(i,:), K, type, p, noimat(i,:), 'bisec');
	%end

	% deconvolved
	w = wienforwd(testymat(i,:), K, type, p, noimat(i,:), sc_a, rho, method);
	decf(i,:) = iwtrans(w, type, p);

	decw(i,:) = decwien(testymat(i,:), K, noimat(i,:), 1);
	deca(i,:) = decall(testymat(i,:), K, 'allp', 1);
	
	% sum of those
	testysum = testysum + testymat(i,:);
	noisum = noisum + noimat(i,:);
	decfsum = decfsum + decf(i,:);
	decwsum = decwsum + decw(i,:);
	decasum = decasum + deca(i,:);
end



	

