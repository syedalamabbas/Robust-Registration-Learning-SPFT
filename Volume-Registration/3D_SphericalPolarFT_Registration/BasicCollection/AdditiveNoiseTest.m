
clear
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


theta_z1 = 2.675*pi/8; % in degrees = 40.5000  >> 0
theta_y = 1.823*pi/8; % in degrees = 63  >> 0
theta_z = 1.84*pi/8; % in degrees = 86.4000  >> 0

surfaceVal = 50;
degree = 31;
isPlotting = 1;
% scale = .3;
scale = .7;
origin=[0 0 0];
vxsize =[1 1 1];

%% Load the original volume data
InDataDir = folderName;

% Without noise
DataFolder = 'benchmark\3DObject_Dataset\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];

inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};

fileIndex = 12;             % Fixed as a dog
SNRdBs = 20:5:40;
F = length(SNRdBs);

all_alphas_    = zeros(1,F);
all_betas_     = zeros(1,F);
all_gammas_    = zeros(1,F);

viewangle1 = -243;
viewangle2 = 23;

color1 = [1 .1 .1];
edgecolor1 = [1 .6 .6]; % dark red
color2 = [.1 1 .1];
edgecolor2 = [.6 1 .6]; % dark green
isLighting_Labels = 1;
voxelSize = 75;
isPlotting = 1;

for ii=1:F   % This loop goes through all the noise levels and computes the rotation errror and RMSE for a fixed object and rotation
    SNRdB = SNRdBs(ii);  % Current noise level in the vertex positioning
    %% Loading the Princeton Database Specific Object
    fileName = [InDataDir   DataFolder '\' inFiles(fileIndex).name];
    loadedFile = load(fileName);
    vertices_noiseless = loadedFile.vertices;
    faces = loadedFile.faces;
    
    if(isPlotting)
        %% Before voxelization
        PlotSurface1(vertices_noiseless,faces); 
        view(viewangle1,viewangle2);
    end
    
    %% Adding noise with the given SNRs for all three coordinates
    vertices_noisy = AddAWGN_XYZPoints( vertices_noiseless, SNRdB );
    if(isPlotting)
        %% Before voxelization
        PlotSurface1(vertices_noisy,faces);
        view(viewangle1,viewangle2);
    end
    
   
    %% After Voxelization
    [bim] = VerticesFacesToVolume( vertices_noiseless, loadedFile.faces, voxelSize );
    if(isPlotting)
        [vertices, faces] = gen_surf_data(bim,origin,vxsize);
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
    end
    volumeOriginal = surfaceVal *double(bim);
    
    [bim] = VerticesFacesToVolume( vertices_noisy, loadedFile.faces, voxelSize );
    if(isPlotting)
        [vertices, faces] = gen_surf_data(bim,origin,vxsize);
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
    end
    volumeOriginalNoisy = surfaceVal *double(bim);
    
    %% The transformation
    t = eye(4);
    t(1:3,1:3)= GetFullRotationMatrixZYZ( theta_z1, theta_y, theta_z );
    tform = affine3d(t);
    
    %% Apply the transformation and display
    volumeRotated = imwarp(volumeOriginalNoisy,tform);     % Rotation with pure signal
    [vertices, faces] =  gen_surf_data(logical(volumeRotated/surfaceVal),origin,vxsize);
    
     isPlotting = 1;
    if(isPlotting)
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
    end
    %% Make Volume Sizes Equal
    [ final_volume1, final_volume2, ~ ] = Pad_MakeVolumeDimensionOdd( volumeOriginal, volumeRotated );
    %% Scale it for speed
    ScalingMat = [scale   0     0     0
        0   scale   0     0
        0     0    scale  0
        0     0      0    1];
    tformTemp =  affine3d(ScalingMat);
    volumeOriginal = imwarp(final_volume1 ,tformTemp);
    volumeRotated = imwarp(final_volume2 ,tformTemp);
    [N, ~, ~] = size(volumeOriginal );
    
    if(~mod(N,2))  % Precautionary measure for not odd case
        volumeOriginal = padarray(volumeOriginal, [1 1 1],'pre');
        volumeRotated = padarray(volumeRotated, [1 1 1],'pre');
    end
    [N, ~, ~] = size(volumeOriginal);  % Now definitely N = odd
    
    B = N+5;
    %% Use the SOFT technique combined with Spherical Polar Fourier Transform
    OutputMatrix = ComputeSOFTRotation_SphericalPolarFT( volumeOriginal, volumeRotated, theta_z1, theta_y, theta_z, B, degree, isPlotting );
    
    %% Find the RMSE for matched Spectral Radius
%     [ final_alpha, final_beta, final_gamma, final_RMSE ] = FindLowestRMSE_Estimate( volumeOriginal, volumeRotated, OutputMatrix, isPlotting );
    final_alpha = OutputMatrix(2,1);   % Picking on the 2nd radius 
    final_beta  = OutputMatrix(2,2);
    final_gamma = OutputMatrix(2,3);
%     %% Show Final Transformation
%     isPlotting = 1;
%     t = eye(4);
%     t(1:3,1:3)= GetFullRotationMatrixZYZ( final_alpha, final_beta, final_gamma );
%     tform = affine3d(t);
%     volumeRotated_GivenTransformation = imwarp(volumeOriginal,tform);     % Rotation with pure signal
%     [ volumeRotated_GivenTransformation, volumeRotated, ~ ] = Pad_MakeVolumeDimensionOdd( volumeRotated_GivenTransformation, volumeRotated );
%     if(isPlotting)
%         %% Plot and Visualize- One way
%         %         [vertices, faces] =  gen_surf_data(logical(volumeRotated/surfaceVal),origin,vxsize);
%         %         PlotSurface1(vertices,faces);
%         %         view(viewangle1,viewangle2);
%         %
%         %         [vertices, faces] =  gen_surf_data(logical(volumeRotated_GivenTransformation/surfaceVal),origin,vxsize);
%         %         PlotSurface1(vertices,faces);
%         %         view(viewangle1,viewangle2);
%         
%         
%         figure,
%         hold on
%         [vertices, faces] =  gen_surf_data(logical(volumeRotated_GivenTransformation),origin,vxsize);
%         CustomPlotVoxelSurface( vertices,faces,  color2, edgecolor2, 0);
%         [vertices, faces] =  gen_surf_data(logical(volumeRotated),origin,vxsize);
%         CustomPlotVoxelSurface( vertices,faces,  color1, edgecolor1, isLighting_Labels  );
%         hold offComputeSOFTRotation_SphericalPolarFT
%         title('Combined directly')
%         
%         
%         
%     end
    
    all_alphas_(ii)   = final_alpha;
    all_betas_(ii)    = final_beta;
    all_gammas_(ii)   = final_gamma;
        
    disp(['Finished Processing, File #', num2str(ii)]);
    
end

all_alphas_

all_betas_

all_gammas_


%% Output of the above program run once

% all_alphas_ =
% 
%     1.9050    2.0053    1.0361    2.0721    1.0026    2.0721    1.0361
% 
% all_alphas_ = [     (pi- 1.9050)    (pi-2.0053)    1.0361    (pi-2.0721)    1.0026    (pi-2.0721)    1.0361];
% all_betas_ =
% 
%     0.6266    0.6601    2.4314    0.7102    2.4314    0.7102    2.4314
% 
% all_betas_ = [   0.6266    0.6601    (pi -2.4314)    0.7102    (pi-2.4314)    0.7102    (pi-2.4314) ]
% all_gammas_ =
% 
%     2.6737    2.5400    5.6148    2.4397    5.5813    2.4397    5.6148
% 
% all_gammas_ = [   (pi- 2.6737)    (pi -2.5400)    (2*pi-5.6148)    (pi-2.4397)    (2*pi-5.5813)    (pi-2.4397)    (2*pi-5.6148) ]
%
% abs(all_alphas_ - theta_z1)
% abs(all_betas_ - theta_y)
% abs(all_gammas_ - theta_z)
% sqrt((1/3)*((abs(all_alphas_ - theta_z1)).^2+(abs(all_betas_ - theta_y)).^2+(abs(all_gammas_ - theta_z)).^2))