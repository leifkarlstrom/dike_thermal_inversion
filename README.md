# dike_thermal_inversion

Scripts associated with forward and inverse modeling of thermal conduction in 1d around a planar dike. 

Inversion is based on 'MCMCHammer', a markov chain monte carlo code written by Kyle Anderson: https://gitlab.com/kander/mcmc_hammer

The forward model is standalone matlab code, finite differences on non-uniform grid with multicomponent (dike+host rock) conduction and partial melting. 

This model has been applied to inverting the following datasets: (U-Th)/He partial resetting in Ap and Zr, Bt/Ar partial resetting, paleomagnetic resetting, Ap fission track. Details are found in 

Goughnour R. L., K. E. Murray, L. Karlstrom, J. Biasi, S. E. Cox, P. O’Sullivan , and B. Finney (2025) "Sensitivity of thermal histories to diverse data constraints: co-inversion of multiple geochronometers and paleomagnetic thermometer at a Columbia River Basalt feeder dike", in press at Geophere.



