function this = restF(f,subset)
f
ff = fft(f);
this = ff(subset);