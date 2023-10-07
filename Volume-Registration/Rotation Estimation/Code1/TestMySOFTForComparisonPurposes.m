clear
clc
close all

%% Load the original volume data
InDataDir = '..\';

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
%% Create holes in the volume
N = 55;
for s = 1:N
    J = bim(:,:,s);
    randomIntegerL = randi(24);
    randomIntegerR = randi(24);
    J((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR,(N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR) = zeros(length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR),length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR));
    bim(:,:,s) = J;
end


[vertices, faces] =  gen_surf_data(bim,origin,vxsize);

PlotSurface1(vertices,faces);
view(220,12);
volumeOriginalNoisy = surfaceVal *double(bim);
%% Add some paths
addpath('..\..\3DImageRegistration');
addpath(genpath('..\..\RadonProject\NewParallelWork'));
addpath(genpath('..\..\MEX_SOFT_Project'));

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

%% Case 2 
volumeOriginal = padarray(volumeOriginal, [19 19 19],'both');
size(volumeOriginal)

volumeRotated = padarray(volumeRotated , [8 6 1],'both');
% volumeRotated = padarray(volumeRotated , [1 0 0],'pre');
size(volumeRotated )

 

 
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
N = 63;
 
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

ComputeRotationSOFTVolumes( volumeOriginal, volumeRotated, theta_z1, theta_y, theta_z, B , degree);

% %% After Computing rotations upto radius = 10 shortlist the number of figures generated and save in following order
% 
% print -dpdf SphericalImageOriginal_1
% print -dpdf SphericalImageNoisyTrans_1
% print -dpdf NoisyTransReg_1
% print -dpdf NoisyTransReg_2
% print -dpdf NoisyTransReg_3
% print -dpdf NoisyTransReg_4
% print -dpdf NoisyTransReg_5
% print -dpdf NoisyTransReg_6
% print -dpdf NoisyTransReg_7
% print -dpdf NoisyTransReg_8
% print -dpdf NoisyTransReg_9
% print -dpdf NoisyTransReg_10