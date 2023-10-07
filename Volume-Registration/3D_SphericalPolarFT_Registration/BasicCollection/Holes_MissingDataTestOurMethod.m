clear all
clc
close all

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
addpath(genpath('..\..\..\RadonProject\NewParallelWork\MEXArrayFireCUDA'))

%% Adding path to Special Correlation On Sphere using SOFT project - My modified form of original and external code
addpath(genpath('..\..\..\MEX_SOFT_Project'))


%% Load the original volume data
InDataDir = folderName;

% Without noise
% DataFolder = 'benchmark\3DObject_Dataset\';
% Templatenamedir=[InDataDir DataFolder '\Template\' ];

% With noise
DataFolder = 'benchmark\3DObject_Dataset\3D_Object_Noise\'; 
Templatenamedir=[InDataDir DataFolder '\TemplateNoise\' ];

inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
for ii=1:length(inFiles)
    inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
end
i = 1;
voxelSize = 75;  % Keep this fixed 

load(inNames_obj{i});
[pa,name,ex]=fileparts(inNames_obj{i}); 

% Before voxelization
PlotSurface1(vertices,faces);
view(220,12);
 
% After Voxelization
% [bim] = verticestovolumefunc(vertices,faces);
[bim] = VerticesFacesToVolume(vertices,faces,voxelSize);
origin=[0 0 0];
vxsize =[1 1 1];
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);
PlotSurface1(vertices,faces);
view(220,12);
surfaceVal = 50;
volumeOriginal = surfaceVal *double(bim);
%% Create holes in the volume

% %% > 80 %
% maxAddZerosStage =1;
% offsetOfVoxels = 23;
% %% > 60
% maxAddZerosStage =2;
% offsetOfVoxels = 20;
% %% > 40
% maxAddZerosStage =3;
% offsetOfVoxels = 17;
% %% > 30
% maxAddZerosStage =3;
% offsetOfVoxels = 15;
%% > 20
maxAddZerosStage =4;
offsetOfVoxels = 12;


[ bim, OverlapPercentage ] = AddHolesInVolume( maxAddZerosStage, offsetOfVoxels, bim );
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);

PlotSurface1(vertices,faces);
view(220,12);
volumeOriginalNoisy = surfaceVal *double(bim);

%% Case 1
% theta_z1 = 1.8*pi/8; % in degrees = 40.5000  >> 0
% theta_y = 2.8*pi/8; % in degrees = 63  >> 0
% theta_z = 3.84*pi/8; % in degrees = 86.4000  >> 0

%% Case 2
theta_z1 = 5.8*pi/8; % in degrees = 40.5000  >> 0
theta_y = 4.8*pi/8; % in degrees = 63  >> 0
theta_z = 3.84*pi/8; % in degrees = 86.4000  >> 0

%% The transformation
t = eye(4);
t(1:3,1:3)= GetFullRotationMatrixZYZ( theta_z1, theta_y, theta_z );
tform = affine3d(t);


%% Apply the transformation and make equal sizes volumes
volumeRotated = imwarp(volumeOriginalNoisy,tform);     % Rotation with pure signal

[vertices, faces] =  gen_surf_data(logical(volumeRotated/surfaceVal),origin,vxsize);

PlotSurface1(vertices,faces);
view(-47,-30);
 

% %% Case 1
% volumeOriginal = padarray(volumeOriginal, [20 20 20],'both');
% size(volumeOriginal)
% 
% volumeRotated = padarray(volumeRotated , [5 7 0],'both');
% volumeRotated = padarray(volumeRotated , [1 0 0],'pre');
% size(volumeRotated )
% 

% %% Case 2 
% volumeOriginal = padarray(volumeOriginal, [19 19 19],'both');
% size(volumeOriginal)
% 
% volumeRotated = padarray(volumeRotated , [8 6 1],'both');
% % volumeRotated = padarray(volumeRotated , [1 0 0],'pre');
% size(volumeRotated )

 
 %% Make Volume Sizes Equal
    [ volumeOriginal, volumeRotated, N ] = Pad_MakeVolumeDimensionOdd( volumeOriginal, volumeRotated );
 
%% Final scaling 
% scale = .66;  % Case 1
scale = .67;
ScalingMat = [scale   0     0     0
                0   scale   0     0
                0     0    scale  0
                0     0      0    1];
tformTemp =  affine3d(ScalingMat);
volumeOriginal = imwarp(volumeOriginal ,tformTemp);
volumeRotated = imwarp(volumeRotated ,tformTemp);
%  N = 63;
  %% Make Volume Sizes Equal
    [ volumeOriginal, volumeRotated, N ] = Pad_MakeVolumeDimensionOdd( volumeOriginal, volumeRotated );
%% Use the SOFT technique combined with Spherical Polar Fourier Transform
B = N+5;
degree = 31;
% % why do we need these two lines ?? 
% volumeOriginal = permute (volumeOriginal , [2 1 3]);   
%  volumeRotated = permute (volumeRotated , [2 1 3]);    

%  %% Save PDF files ?
% % print -dpdf Original3DSurfaceData  -- Figure 1 select
% % print -dpdf Voxelized3DSurfaceData  -- Figure 2 select
% % print -dpdf NoisyVoxelized3DSurfaceData  -- Figure 3 select
% % print -dpdf NoisyTransformedVoxelized3DSurfaceData  -- Figure 4 select
isPlotting = 1;
OutputMatrix = ComputeSOFTRotation_SphericalPolarFT( volumeOriginal, volumeRotated, theta_z1, theta_y, theta_z, B, degree, isPlotting );


%% Output of the above script for different overlaps --- Our Method
% > 80% 
% R_actual  = [2.2777    1.8850    1.5080];
% R_Estimated = [(pi-0.8378)    1.8762    1.6755];  % At Spectral Radius 10
% rms(R_actual-R_Estimated) % RMSE computed 
% sqrt(sum((R_actual-R_Estimated).^2)/3)  = 
% > 60% 
% R_actual  = [2.2777    1.8850    1.5080];
% R_Estimated = [(pi-0.8378)   1.8762     1.695];  % At Spectral Radius 10
% rms(R_actual-R_Estimated) % RMSE computed 
% > 40% 
% R_actual  = [2.2777    1.8850    1.5080];
% R_Estimated = [(pi-0.8378)    1.8937    1.6755];  % At Spectral Radius 10
% rms(R_actual-R_Estimated) % RMSE computed 
% > 30% 
% R_actual  = [2.2777    1.8850    1.5080];
% R_Estimated = [(pi-0.8029)    1.8762    1.6755];  % At Spectral Radius 10
% rms(R_actual-R_Estimated) % RMSE computed 
% > 20% 
% R_actual  = [2.2777    1.8850    1.5080];
% R_Estimated = [(pi-0.8727)    1.9286    1.1868];  % At Spectral Radius 10
% rms(R_actual-R_Estimated) % RMSE computed 




%% Output of the Method used in TIP 2013 paper - 
% > 80% 
% rms(R_actual-R_Estimated) % RMSE computed 
 