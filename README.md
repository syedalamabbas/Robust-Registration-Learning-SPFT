# Robust-Registration-Learning-SPFT
The source code for robust registration and learning using multi-radii spherical polar Fourier transform

General comment: most of the folders should have readme files that describe the contents and purpose of it.

Overall it is:

1. SPFT-ArrayFire has C++ and ArrayFire library code to simply compute the SPFT given volumes.
2. Inside Machine Learning folder it is all organized as follows:
      i. ModelNet40 preprocessing contains the scripts necessary to obtain just the magnitude spectrum from multi-radii SPFT. 
      ii. It also contains a helpful script to create tensorflow record that is used to train and test the classification network.
      iii. Common includes some useful and helpful MATLAB scripts extracted from the conventional volumetric registration algorithm.
3. For Volume Registration:
      i. SOFT2.0 contains the base C code for computing spherical correlations given two equiangular grids
      ii. SOFT_MEX is a project that allows for create of MEX file to be used in MATLAB with the registration scripts.



Acknowledgements(each of the software or datasets should likely have their respective licenses in its folder):

1. SOFT: SO(3) Fourier Transforms
2. ModelNet40 dataset
3. 
