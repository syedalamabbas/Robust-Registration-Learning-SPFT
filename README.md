# Robust-Registration-Learning-SPFT
#+TITLE: Robust Registration and Learning using Multi-radii Spherical Polar Fourier Transform.

[[file:SPFTMotivationCat.png]
[[file:MotivationCatSPFT.png]

General comment: Most of the folders should have readme files that describe the contents and purpose of it.

Overall it is all organized as follows:

1. SPFT-ArrayFire has C++ and ArrayFire library code, which is common for both volumetric registration as well as for machine learning, to simply compute the SPFT given volumes. However it is only for the machine learning work, involving tensorflow, python and jupyter   that we had to extract this. Much of the experiments for volumetric registration were conducted in MATLAB where both SPFT and SOFT C/C++ projects had to be interfaced as MEX libraries. 
2. For Machine Learning:
      i. ModelNet40 preprocessing contains the scripts necessary to obtain just the magnitude spectrum from multi-radii SPFT. 
      ii. It also contains a helpful script to create tensorflow record that is necessary and used to train and test the classification  network, the slightly modified spherical CNN. Much of the data examples, both in archive as well as the extracted images, had to be erased for the repository, otherwise it was bloating it over many GBs and causing git push unexpected hung-up failures. 
      iii. Common includes some useful and helpful MATLAB scripts extracted from the conventional volumetric registration algorithm for separately visualizing the images in this setting.
3. For Volume Registration:
      i. SOFT2.0 contains the base C code for computing spherical correlations given two equiangular grids on a single sphere.
      ii. SOFT_MEX is a project that allows for create of MEX file to be used in MATLAB with the registration scripts, this helps in visualizing the correlation results nicely in MATLAB.



Acknowledgements, in no particular order (each of the software or datasets should likely have their respective licenses in its folder):

1. SOFT: SO(3) Fourier Transforms
2. ModelNet40 dataset
3. 
