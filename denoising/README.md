# De-noising images 

4 local and 1 nonlocal denoising algorithms

# Requirements

MATLAB with image processing toolbox

OR,

OCTAVE with image package for octave like this:

0. Install octave
1. Download the image package from https://octave.sourceforge.io/packages.php
2. Run octave with `octave`.
3. Install the package `pkg install /path/to/image-pkg.tar.gz`
4. Load package with `pkg load image`

# Usage

If you don't have a noisy image, create on by `imnoise()`
* Gaussian denoising (blur): `gaussF2(nbdsize, variance, image_handle)`
* Median denoising (blur): `medF(nbdsize, imgage_handle)`
* Susan filter: `susanF2(nbdsize, variance, image_handle)`
* Yaroslavsky Filter: `ynf(nbdsize,variance, image_handle)`

* Nonlocal Mean: `NL3(nbdsize, h, variance, image_handle)` (Only nonlocal algorithm)

Note: To apply these methods to a color image, you need to apply it to each channel (out of R,G,B) individually and merge them in the end. See `handleCol.m` for example.

# Reference

* "Deblurring and Denoising of Images by Nonlocal Functionals“ by Stefan Kin-dermann, Stanley Osher, and Peter W. Jones.
* "On image denoising methods" by A. Buades, B. Coll, and J. M. Morel,  preprint,CMLA, Cachan, France, 2004.
