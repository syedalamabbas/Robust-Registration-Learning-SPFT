clear all;
close all;
clc;

%% Finding the translation of the 3D models

%% Adding path to Spherical Harmonics External code
% - Many of the function in it are also used here
folderName = '..\..\..\Rotation Estimation\';
addpath(genpath(folderName));


%% Adding path to the root folder of this project i.e. 3D_SphericalPolarFT_Registration
addpath(genpath('..\'))

%% Adding path to Spherical Polar Fourier Transform 2016 TSP - My own external code
addpath(genpath('..\..\..\RadonProject\NewParallelWork\NVIDIA_3DSphericalDFT'))
addpath(genpath('..\..\..\RadonProject\NewParallelWork\NVIDIA_2DPolarDFT'))
addpath(genpath('..\..\..\RadonProject\NewParallelWork\NVIDIA_1DFrFT'))

%% Adding path to Special Correlation On Sphere using SOFT project - My modified form of original and external code
addpath(genpath('..\..\..\MEX_SOFT_Project'))

%% Load the original volume data
InDataDir = folderName;

% Without noise
DataFolder = 'benchmark\3DObject_Dataset\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];

% % With noise
% DataFolder = 'benchmark\3DObject_Dataset\3D_Object_Noise\';
% Templatenamedir=[InDataDir DataFolder '\TemplateNoise\' ];

inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
for ii=1:length(inFiles)
    inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
end
i = 1;

load(inNames_obj{i});
[pa,name,ex]=fileparts(inNames_obj{i});

% Before voxelization
PlotSurface1(vertices,faces);
view(220,12);

% After Voxelization
[bim] = verticestovolumefunc(vertices,faces);
origin=[0 0 0];
vxsize =[1 1 1];
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);
PlotSurface1(vertices,faces);
view(220,12);
surfaceVal = 50;
volumeOriginal = surfaceVal *double(bim);
% %% Create holes in the volume
% N = 55;
% % for s = 1:N
% %     J = bim(:,:,s);
% %     randomIntegerL = randi(24);
% %     randomIntegerR = randi(24);
% %     J((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR,(N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR) = zeros(length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR),length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR));
% %     bim(:,:,s) = J;
% % end
%
% [vertices, faces] =  gen_surf_data(bim,origin,vxsize);
% PlotSurface1(vertices,faces);
% view(220,12);
% volumeOriginalNoisy = surfaceVal *double(bim);

%% The transformation and application
t = eye(4);
translation_actual = [58 43 26];
alpha = 2.345;    % This solution is rotationally invariant ???
beta  = 1.234;
gamma = 1.7645;
t(1:3,1:3) = eye(3);GetFullRotationMatrixZYZ(alpha, beta, gamma);
t(4,1:3)= translation_actual;
tform = affine3d( t);
RA = imref3d(size(volumeOriginal));

[volumeTranslated, RB] = imwarp(volumeOriginal,RA, tform);     % Translation with pure signal
% volumeTranslated = padarray(volumeOriginal, translation_actual, 'pre');


%% Apply the transformation and make equal sizes volumes
[vertices, faces] =  gen_surf_data(logical(volumeTranslated/surfaceVal),origin,vxsize);
PlotSurface1(vertices,faces);
view(220,12);
 
%% Make Volume Sizes Equal
[ volumeOriginal, volumeTranslated, N ] = Pad_MakeVolumeDimensionOdd( volumeOriginal, volumeTranslated );

%% Computing the translation 
[ translation_vector ] = ComputeTranslationWithVolumes( volumeOriginal, volumeTranslated );

%% Compare
Error = translation_actual - translation_vector

