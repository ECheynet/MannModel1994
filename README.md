# MannModel1994
Numerical implementation of the uniform-shear model by Mann (1994)

[![View Uniform shear model (Mann, 1994) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/67055-uniform-shear-model-mann-1994)

The computation of the one-point spectra, co-spectra and coherence of turbulent wind is conducted using the uniform-shear model from Mann (1994) [1]. The goal is to describe the spatial structure of stationary homogeneous turbulence under a neutral atmospheric stratification using only 3 adjustable parameters.
The present submission contains:
- The function MannTurb.m that computes the sheared spectral tensor
- The function MannCoherence that computes the wind co-coherence
- A LiveScript for the example file
- A data file GreatBeltSpectra.mat that contains the wind spectra from the Great belt bridge experiment for comparison with the computed spectra.

This is the first version of the submission. Some bugs may still exist. Any question, comment or suggestion is warmly welcomed.

References

[1] Mann, J. (1994). The spatial structure of neutral atmospheric surface-layer turbulence. Journal of fluid mechanics, 273, 141-168.
