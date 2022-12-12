# Numerical implementation of the uniform-shear model
Matlab implementation of the uniform shear model by Mann (1994)

[![View Uniform shear model (Mann, 1994) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/67055-uniform-shear-model-mann-1994)
[![DOI](https://zenodo.org/badge/249148606.svg)](https://zenodo.org/badge/latestdoi/249148606)
[![Donation](https://camo.githubusercontent.com/a37ab2f2f19af23730565736fb8621eea275aad02f649c8f96959f78388edf45/68747470733a2f2f77617265686f7573652d63616d6f2e636d68312e707366686f737465642e6f72672f316339333962613132323739393662383762623033636630323963313438323165616239616439312f3638373437343730373333613266326636393664363732653733363836393635366336343733326536393666326636323631363436373635326634343666366536313734363532643432373537393235333233303664363532353332333036313235333233303633366636363636363536353264373936353663366336663737363737323635363536653265373337363637)](https://www.buymeacoffee.com/echeynet)

## Summary
The computation of the one-point spectra, co-spectra and coherence of turbulent wind is conducted using the uniform-shear model from Mann (1994) [1]. The goal is to describe the spatial structure of stationary homogeneous turbulence under a neutral atmospheric stratification using only 3 adjustable parameters.

## Content

The present submission contains:
- The function MannTurb.m that computes the sheared spectral tensor
- The function MannCoherence that computes the wind co-coherence
- A LiveScript for the example file
- A data file GreatBeltSpectra.mat that contains the wind spectra from the Great belt bridge experiment for comparison with the computed spectra.

Any question, comment or suggestion is warmly welcomed.

## References

[1] Mann, J. (1994). The spatial structure of neutral atmospheric surface-layer turbulence. Journal of fluid mechanics, 273, 141-168.
