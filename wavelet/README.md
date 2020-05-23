# Wavelet

Wavelet transforms and related tools in 1 dimension

# Requirements

OCTAVE or MATLAB

# Usage

Main functions:
* `[u, v] = filt(type, p)` returns the parent wavelets `u`, `v` used in wavelet transform used of type `type`.
* `w = wtrans(z, type, p)` returns the `p`th stage wavelet transform of the vector `z` using the parent wavelets of type `type`.
* `z = iwtrans(w, type, p)` returns the `p`th stage _inverse_ wavelet transform of the vector `w` using the parent wavelets of type `type`. If the `type` and `p` are the same, `iwtrans` and `wtrans` should be inverses of each other upto some _small_ error.
* `bases(type, p, N)` to plot the bases with respect to (the translations of) which the wavelet transform is computed. The implementation is not recursive.
* `plotbasis(type, p, N)` (better) plots the basis elements by first computing the full basis matrix using `getbaismat`
* `plotcoeffs(w,p)` plots the coefficients of wavelet transform `w`, the `p`-th stage wavelet transform
* `coeffs(w, p, q)` prints the q-th level wavelet coefficients of the p-th resolution wavelet transform. `q` can be from 1 to `p+1`. q=p+1 represents the coarsest wavelet level.

* [obsolete] `isolate(z, type, p)` plots various wavelet levels of the p-th resolution wavelet deconvolution of the signal `z`

Compression and error
* `keeplarge` zeros out smaller values
* `w = compress(z, type, p, K)` returns 
* `relerr(z, typelist, kMax, p, q)` returns the relative error matrix computed in `q`-norm ...
* `compareErr(z, typelist, kMax, p, normlist)`

Parameters:
* `type` can be 'shan' for Shannon's wavelets, 'd`n`' for Daubechies wavelets where `n` can be 2, 4, ... , 20.
* `p` should be a power of 2 (preferably). If `type`= 'shan', `p` must divide the length of `z` or `w`.

Other tools (codes are self-explanatory):
* `filt` contains the parent wavelets of various wavelet bases
* `wrec`, `iwrec` are recursive implementation of wavelet transform and inverse wavelet transforms, also called Fast Wavelet Transform. Works best (fastest) on vectors of length of type 2^n for some natural number n.
* `realconv(a,b)` is a convolution that returns real number
* `fold(z)` folds a vector in half
* `up`, `down` does upsampling and downsampling of a signal by zeros
* `v = getother(u)` given one parent wavelet `u`, returns the other `v` using the conjugate mirror filter lemma
* `B = getbasismat(type, p, N)` generates a matrix of dimension `(p+1)xN` whose j-th row contains the basis element which is translated to generate the vectors used to compute the j-th level wavelet coefficients of a p-th resolution wavelet transform of a signal of length `N`.


Fourier based Deconvolution
* `[fw, mult] = fdecwien(fsig, fimp, noise, scaling)` returns the Fourier transform `fw` of the Wiener deconvolution with the (usually unknown) signal strength as ratio of the supplied signal strength and the variance of the impulse response. `mult` returns the Fourier shrinkage parameter used in the deconvolution.

Wavelet based deconvolution
* `alpha = getoptsc(z, K, type, p, sigma, rootmethod)` compute the optimal  scaling parameter vector `alpha` of lenth `p+1` that minimizes the error in ForWaRD algorithm. Here `rootmethod` can be either `search` (for searching through uniformly located numbers between 0 and 1) or `bisec` (for a bisection method, much faster and accurate)

ANITA related tools
* `getdata` generates some test data to play with, both theoretical and recorded by ANITA
* `plotfanita(z)` plots the absolute value of the Fourier transform of an ANITA signal in the frequency domain in MHz unit
* `plotfanitaS(z)` (with a different sample rate) plots the absolute value of the Fourier transform of an ANITA signal in the frequency domain in MHz unit
* `getantenna(num)` outputs the filename of the impulse response of the antenna, given the absolute index of the antenna
* `[snrval, M, m, a, b] = getsnr(z)` computes the signal-to-noise ratio of a temporally localized signal `z`. Here, `M` and `m` are the max and the min of the signal, the time interval `[a,b]` has 10 times the length of the peak region and is a neighborhood of the peak region.
* `wpxsc = getwpxsc(type, p, sigma)` returns the optimal scaling values `alpha` required for the deconvolving the ANITA signal `wpx` (see `getdata` for the test signal), using the bisection method 

Other tools
* `z = padfreq(w, N, eps)` increases the resolution of the signal (i.e. upsamples) `w` by zero-padding followed by smoothing. If `N > length(w)` it adds zeros in the frequency domain. However, the process can add Gibbs phenomenon in the time domain if the frequency does not decay to zero. To avoid this, we apply a smooth Fourier cut-off function (a Plank-taper window in the frequency domain) with `eps` (between 0 and 1) being the proportion of the _actual_ frequency of the signal `w` to get scaled down. For example, `eps = 1` scales the entire frequency range of the original signal. 
To apply the Plank-taper window without padding, use `N = length(w)` and any value of `eps`.
* `planktaper(N,eps)` generates a Plank-taper window of length `N` and tapering window of length `eps*N`. E.g. `eps = 0.5` the function obtains value 1 on a set of measure zero
* `croscor(f,q)` computes the circular cross-correlation between two vectors of same length
* `cpuow(z)` computes the cumulative power distribution of a signal `z`
* `deriv(z)` computes the difference of consecutive terms (a crude derivative)


# Reference

* M. Frazier, _An Introduction to Wavelets through Linear Algebra_
* Thanks to Dr. Svetlana Roudenko for teaching the wavelets
* Thanks to Dr. Peter Gorham for providing support and ANITA data and Manuel Olmedo for making the data matlab-friendly
