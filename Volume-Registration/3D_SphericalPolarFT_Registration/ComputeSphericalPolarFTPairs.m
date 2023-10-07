function [ VolumeSphericalGrid, VolumeRotatedSphericalGrid ] = ComputeSphericalPolarFTPairs( volumeOriginal,volumeRotated, B )
%COMPUTESPHERICALPOLARFTPAIRS
% This function simply computes the Spherical Polar FT Pairs given two
% volumes
% volumeOriginal             - The first N^3 volume that is the reference
% volumeOriginal             - The second N^3 volume that is the rotated volume w.r.t first volume
%====================================================================
% Written on June 22, 2016 by Syed Alam Abbas.
% Revised on March 12th, 2022 by Syed Alam Abbas.
%====================================================================

%% Computing the Spherical Polar Fourier Transform using my own technique
%%% why do we need these two lines ?? Not sure something to do with SOFT Code rotation conventions
isUsingGPU = 0;
volumeOriginal = permute (volumeOriginal , [2 1 3]);
volumeRotated = permute (volumeRotated , [2 1 3]);

noOfAnglesTheta               = 2*B;    % Must be even  Special Conversion and the ratio
noOfAnglesPhi                 = B;      % Must be even
[~,~,N] = size(volumeOriginal);

if(~isUsingGPU)     % Then use CPU
    tic
    VolumeSphericalGrid        = abs( VectorizedCompute3DSphericalDFT( (volumeOriginal),  noOfAnglesTheta, noOfAnglesPhi));       % Killing phase information grabbing magnitude only
    toc
    tic
    VolumeRotatedSphericalGrid = abs( VectorizedCompute3DSphericalDFT( (volumeRotated),  noOfAnglesTheta, noOfAnglesPhi));        % Killing phase information grabbing magnitude only
    toc
else
    gVol1 = gpuArray(volumeOriginal);
    tic
    [realImage, ImagImage]= SphericalPolar3DTransformMEXArrayFire(gVol1, noOfAnglesTheta, noOfAnglesPhi,N, 0);
    VolumeSphericalGrid = gather(abs(complex(realImage,ImagImage)));
    toc
    clear realImage
    clear ImagImage
    clear gVol1
    clear SphericalPolar3DTransformMEXArrayFire
    g = gpuDevice(1);
    reset(g);
    
    gVol2 = gpuArray(volumeRotated);
    tic
    [realImage, ImagImage]= SphericalPolar3DTransformMEXArrayFire(gVol2, noOfAnglesTheta, noOfAnglesPhi,N, 0);
    VolumeRotatedSphericalGrid = gather(abs(complex(realImage,ImagImage)));
    toc
    clear realImage
    clear ImagImage
    clear gVol2
    clear SphericalPolar3DTransformMEXArrayFire
    g = gpuDevice(1);
    reset(g);
    %     'C:\Program Files\MATLAB\R2013b\toolbox\matlab\winfun\winqueryreg.mexw64'
end