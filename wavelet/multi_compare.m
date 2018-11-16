% compare the relative error in different methods

% get data
multi_data;

% get the optimum scaling with fone()
alpha = multi_getsc_fone(wax, aximp, testyori, 'meyer', 5, noiseax);

% use this scaling to compute the error
[z, unthr] = multi_fone(wax, aximp, testyori, 'meyer', 5, noiseax, alpha, 1, 'hard');

fone_err = rele(z, testyori)
