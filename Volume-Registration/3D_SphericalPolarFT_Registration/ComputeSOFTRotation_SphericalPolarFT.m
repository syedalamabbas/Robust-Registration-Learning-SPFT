function OutputMatrix = ComputeSOFTRotation_SphericalPolarFT( volumeOriginal, volumeRotated, theta_z1, theta_y, theta_z, B , degree, isPlotting)
%% ComputeSOFTRotation_SphericalPolarFT calculates the Spherical Polar Fourier
% transform then uses SOFT technique to do correlation on spheres at each
% spectral radius to find the three angles, it is multi-layered and
% exceptionally robust to noise, distortions, blurring etc.
% volumeOriginal             - The first N^3 volume that is the reference 
% volumeOriginal             - The second N^3 volume that is the rotated volume w.r.t first volume
% theta_z1, theta_y, theta_z - These angles are the true transformations only needed for plotting purposes
% B                          - Bandwidth for the SOFT technique
% degree                     - Spherical Harmonics coefficient degrees for SOFT technique
% isPlotting                 - Plotting the solution = 1 or not = 0
% OutputMatrix               - Spectral Radius vs [alpha, beta, gamma]
%====================================================================
% Written on June 8th, 2016 by Syed Alam Abbas.
% Revised on March 12th, 2022 by Syed Alam Abbas.
%====================================================================

%% Computing the Spherical Polar Fourier Transform using my own technique
[ VolumeSphericalGrid, VolumeRotatedSphericalGrid ] = ComputeSphericalPolarFTPairs( volumeOriginal,volumeRotated, B );

%% Computing the Rotation estimation based on SOFT technique for each spectral layer
[ OutputMatrix ] = ComputeSOFTRotationPairs( VolumeSphericalGrid, VolumeRotatedSphericalGrid, theta_z1, theta_y, theta_z, B , degree, isPlotting );
end

