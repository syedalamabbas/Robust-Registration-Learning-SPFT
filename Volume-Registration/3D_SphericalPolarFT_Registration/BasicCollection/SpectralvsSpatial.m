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
% With noise
DataFolder = 'benchmark\3DObject_Dataset\3D_Object_Noise\'; 
Templatenamedir=[InDataDir DataFolder '\TemplateNoise\' ];

inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
for ii=1:length(inFiles)
    inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
    [pa,name,ex]=fileparts(inNames_obj{ii}); 
    ii, name
end

obj_id = 10;
voxelSize = 145;  % Keep this fixed 

'Voxel size='
voxelSize

load(inNames_obj{obj_id});
[pa,name,ex]=fileparts(inNames_obj{obj_id}); 

%Normalizing the input vertices
[r, c] = size(vertices)
avgx = sum(vertices(:, 1))/r;
avgy = sum(vertices(:, 2))/r;
avgz = sum(vertices(:, 3))/r;
vertices(:, 1) = vertices(:, 1) - avgx;
vertices(:, 2) = vertices(:, 2) - avgy;
vertices(:, 3) = vertices(:, 3) - avgz;

% Rotate the object
x = 0; 
y = 0;
z = pi;
Rx = [      1       0       0; ...
            0  cos(x) -sin(x); ...
            0  sin(x)  cos(x)];

Ry = [ cos(y)       0  sin(y); ...
            0       1       0; ...
      -sin(y)       0  cos(y)];

Rz = [ cos(z) -sin(z)       0; ...
       sin(z)  cos(z)       0; ...
            0      0        1];

R = Rz*Ry*Rx;
vertices = vertices * R;

% Before voxelization
PlotSurface1(vertices,faces);
view(146, 5);
 
% After Voxelization
[bim] = VerticesFacesToVolume(vertices,faces,voxelSize);
origin=[0 0 0];
vxsize =[1 1 1];
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);
PlotSurface1(vertices,faces);
view(146, 5);

surfaceVal = 50;
volumeOriginal = surfaceVal *double(bim);

maxAddZerosStage =2;
offsetOfVoxels = 30;

[ bim, OverlapPercentage ] = AddHolesInVolume( maxAddZerosStage, offsetOfVoxels, bim );
[vertices, faces] =  gen_surf_data(bim,origin,vxsize);

PlotSurface1(vertices,faces);
view(146, 5);
volumeOriginalNoisy = surfaceVal *double(bim);

%% Test sphere warping
% [I, map] = imread('forest.tif');
% [X, Y, Z]= sphere(40);
% figure, 
% warp(X, Y, Z, I, map)
% surface(x,y,z,'FaceColor', 'none','EdgeColor',[.7 .7 .7]);
% axis equal

% ImageOnSphere = rgb2gray(imread('out8x8.jpg')) ;

ImageOnSphere = rgb2gray(imread('HP+Syn.jpg'));

figure, 
imagesc(ImageOnSphere);
% colorbar
axis off
% hold off

map = parula(256);
colorImg = ind2rgb(ImageOnSphere, map);
[x,y,z]= sphere(40);
figure,
hold on
warp(x,y,z, colorImg, map );
% surface(x,y,z,'FaceColor', 'none','EdgeColor',[0.7 0.7 0.7]);
hold off
axis equal;
grid on
xlabel('x')
ylabel('y')
zlabel('z')
%% Case 2
theta_z1 = 1.8*pi/8; 
theta_y = 4.8*pi/8; 
theta_z = 3.84*pi/8; 

%% The transformation
t = eye(4);
t(1:3,1:3)= GetFullRotationMatrixZYZ( theta_z1, theta_y, theta_z );
tform = affine3d(t);

volumeRotatedNoisy = imwarp(volumeOriginalNoisy,tform);     % Rotation with pure signal
[vertices, faces] =  gen_surf_data(logical(volumeRotatedNoisy/surfaceVal),origin,vxsize);
PlotSurface1(vertices,faces);
view(121, 66);

%% Test csv input and output
N = 65; % Output is N x N x N volume
test_voxels = randn(N, N, N);
csvwrite('SerializedVolume.csv' ,test_voxels(:));
B = (N -1)/2;
noOfAnglesTheta               = 2*B;  
noOfAnglesPhi                 = B;  
test_sphericalGrid            = abs( VectorizedCompute3DSphericalDFT( test_voxels,  noOfAnglesTheta, noOfAnglesPhi));
csvwrite('SerializedSphericalGrid.csv' ,test_sphericalGrid(:));

%% Plotting spectral signatures
B=72
noOfAnglesTheta               = 2*B;    % Must be even  Special Conversion and the ratio
noOfAnglesPhi                 = B;      % Must be even
[~,~,N] = size(volumeOriginal);

tic
VolumeSphericalGrid             = abs( VectorizedCompute3DSphericalDFT( (volumeOriginal),  noOfAnglesTheta, noOfAnglesPhi));
toc
tic
VolumeSphericalGridNoisy        = abs( VectorizedCompute3DSphericalDFT( (volumeOriginalNoisy),  noOfAnglesTheta, noOfAnglesPhi));  
toc
% tic
% VolumeSphericalGridNoisyRotated = abs( VectorizedCompute3DSphericalDFT( (volumeRotatedNoisy),  noOfAnglesTheta, noOfAnglesPhi)); 
% toc

dcVal = VolumeSphericalGrid(1,1,(N-1)/2+1);
dcValNoisy = VolumeSphericalGridNoisy(1,1,(N-1)/2+1);

for k =1:3
    GetSphericalImage(  VolumeSphericalGrid/dcVal, k , 1);
    GetSphericalImage(  VolumeSphericalGridNoisy/dcValNoisy, k , 1);
%     GetSphericalImage(  VolumeSphericalGridNoisyRotated/dcValNoisy, k , 1);
end

