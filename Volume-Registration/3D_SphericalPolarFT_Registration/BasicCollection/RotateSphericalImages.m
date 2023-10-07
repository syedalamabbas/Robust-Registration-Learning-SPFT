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

obj_id = 10; %Is a cat
voxelSize = 145;  % Keep this fixed

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

    x = 0;
    y = 0;
    
    Rx = [      1       0       0; ...
        0  cos(x) -sin(x); ...
        0  sin(x)  cos(x)];
    
    Ry = [ cos(y)       0  sin(y); ...
        0       1       0; ...
        -sin(y)       0  cos(y)];

K = 20;
vertices_new = vertices;
for k =1:K
    % Rotate object only along z axis
    z = (k-1) * 2 * pi/(K - 1);
    Rz = [ cos(z) -sin(z)       0; ...
       sin(z)  cos(z)       0; ...
            0      0        1];
    vertices_new = vertices * Rz; %Clockwise
    
    % Before voxelization
    PlotSurface1(vertices_new,faces);
    pause(0.001)
end
