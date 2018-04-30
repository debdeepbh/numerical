# Wavelet

Wavelet transforms and related tools in 1 dimension

# Requirements

OCTAVE or MATLAB

# Usage

Main functions:
* `[u, v] = filt(type, p)` returns the parent wavelets `u`, `v` used in wavelet transform used of type `type`.
* `w = wtrans(z, type, p)` returns the `p`th stage wavelet transform of the vector `z` using the parent wavelets of type `type`.
* `z = iwtrans(w, type, p)` returns the `p`th stage _inverse_ wavelet transform of the vector `w` using the parent wavelets of type `type`. If the `type` and `p` are the same, `iwtrans` and `wtrans` should be inverses of each other upto some _small_ error.
* `w = compress(z, type, p, K)` returns 
* `bases(type, p, N)` to plot the bases with respect to (the translations of) which the wavelet transform is computed. The implementation is not recursive.
* `relerr(z, typelist, kMax, p, q)` returns the relative error matrix computed in `q`-norm ...
* `compareErr(z, typelist, kMax, p, normlist)`

Parameters:
* `type` can be 'shan' for Shannon's wavelets, 'd`n`' for Daubechies wavelets where `n` can be 2, 4, ... , 20.
* `p` should be a power of 2 (preferably). If `type`= 'shan', `p` must divide the length of `z` or `w`.

Other tools (codes are self-explanatory):
* `wrec`, `iwrec` are recursive implementation of wavelet transform and inverse wavelet transforms, also called Fast Wavelet Transform. Works best (fastest) on vectors of length of type 2^n for some natural number n.
* `realconv(a,b)` is a convolution that returns real number
* `fold(z)`
* `up`, `down` sampling
* `getother`
* `keeplarge` zeros out smaller values

# Reference

* M. Frazier, _An Introduction to Wavelets through Linear Algebra_
* Thanks to Dr. Svetlana Roudenko for teaching the material
