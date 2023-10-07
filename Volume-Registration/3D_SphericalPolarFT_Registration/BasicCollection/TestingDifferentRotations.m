close all
clear
clc
% Just some additional scripts not used in the paper by Syed Alam Abbas,
% PAMI - Robust $3$D Registration using Spherical Polar Fourier Transform and Spherical Harmonics

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

%% Discretize the angular grids in the three axis rotations
angle_resolution = 10*pi/180;  % 10 degrees
alphas = angle_resolution:angle_resolution: 2*pi-angle_resolution; % In Radians , seperated by 10 degrees each
betas  = angle_resolution:angle_resolution: 2*pi-angle_resolution; % In Radians , seperated by 10 degrees each
gammas = angle_resolution:angle_resolution: 2*pi-angle_resolution; % In Radians , seperated by 10 degrees each

K = length(alphas);

%% Gather the volume on which we operate the estimation procedure
surfaceVal = 50;
fileIndex = 7;
viewangle1 =  220;
viewangle2 = 301;
origin=[0 0 0];
vxsize =[1 1 1];

InDataDir = folderName;
DataFolder = 'benchmark\3DObject_Dataset\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];
inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};

fileName = [InDataDir   DataFolder '\' inFiles(fileIndex).name];
load(fileName);

[pa,name,ex]=fileparts(fileName);
PlotSurface1(vertices,faces);
view(viewangle1,viewangle2);
title('Original')

[bim] = verticestovolumefunc(vertices,faces);
[vertices, faces] = gen_surf_data(bim,origin,vxsize);
PlotSurface1(vertices,faces);
view(viewangle1,viewangle2);
volumeOriginal = surfaceVal *double(bim);
title('Voxelized')

%% Perform the estimation procedure and compute the RMSE for different angles
degree = 35;
isPlotting = 0;

X_Coordinates = zeros(K,K,K,K);
Y_Coordinates = zeros(K,K,K,K);
Z_Coordinates = zeros(K,K,K,K);
Final_RMSEs = zeros(K,K,K,K);
Final_Angular_RMSE = zeros(K,K,K,K);

count =0;
for i=30:K
    for j =9:K
        for k =29:K
            %% Actual rotation based on the grid
            actual_rotation = [alphas(i), betas(j), gammas(k)];
            
            %% Create and Apply the transformation
            t = eye(4);
            t(1:3,1:3)= GetFullRotationMatrixZYZ( actual_rotation(1), actual_rotation(2), actual_rotation(3) );
            tform = affine3d(t);
            volumeRotated = imwarp(volumeOriginal,tform);     % Rotation with pure signal
            [ volumeOriginal, volumeRotated, N ] = Pad_MakeVolumeDimensionOdd( volumeOriginal, volumeRotated );
            if(~mod(N,2))  % Precautionary measure for not odd case
                volumeOriginal = padarray(volumeOriginal, [1 1 1],'pre');
                volumeRotated = padarray(volumeRotated, [1 1 1],'pre');
            end
            [N, ~, ~] = size(volumeOriginal);  % Now definitely N = odd
            B = N+5;
            
            %% Use the SOFT technique combined with Spherical Polar Fourier Transform
            OutputMatrix = ComputeSOFTRotation_SphericalPolarFT( volumeOriginal, volumeRotated, actual_rotation(1), actual_rotation(2), actual_rotation(3), B, degree, isPlotting );
            
            %% Find the RMSE for matched Spectral Radius
            
            [ final_alpha, final_beta, final_gamma, final_RMSE ] = FindLowestRMSE_Estimate( volumeOriginal, volumeRotated, OutputMatrix, isPlotting );
            estimated_rotation = [final_alpha, final_beta, final_gamma];
            
            X_Coordinates(i,j,k) = i;
            Y_Coordinates(i,j,k) = j;
            Z_Coordinates(i,j,k) = k;
            Final_RMSEs(i,j,k) = final_RMSE;
            Final_Angular_RMSE(i,j,k) = rms(actual_rotation - estimated_rotation);
            %             Final_RMSEs(i,j,k) = rand(1,1); % Used Only for testing
            
            PlotTwoRotationAngles( actual_rotation, estimated_rotation );
            count = count + 1;
            print( ['D:\GitHub\Figures\RotationTest', num2str(count) ] ,'-dpdf');
            
        end
    end
end


%% Plotting the final RMSE's in the form of surface plots

% E.g. This is a simple example
% [X,Y] = meshgrid(-8:.5:8);
% R =  sqrt(X.^2 + Y.^2) + eps;
% Z = sin(R)./R;
% mesh(X,Y,Z)

k_index = 24;
X = X_Coordinates(:,:,k_index);
Y = Y_Coordinates(:,:,k_index);
C = Final_RMSEs(:,:,k_index);
figure,
surf(X,Y,C,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
colormap hot
alpha 0.1
title(['RMSE at slice z =',num2str(k_index)])


k_index = 20;
X = X_Coordinates(:,:,k_index);
Y = Y_Coordinates(:,:,k_index);
C = Final_RMSEs(:,:,k_index);
figure,
surf(X,Y,C,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
colormap hot
alpha 0.1
title(['RMSE at slice z =',num2str(k_index)])


