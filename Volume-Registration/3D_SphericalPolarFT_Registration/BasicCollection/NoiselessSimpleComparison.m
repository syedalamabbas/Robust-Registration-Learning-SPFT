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
% addpath(genpath('..\..\..\RadonProject\NewParallelWork\NVIDIA_3DSphericalDFT'))
% addpath(genpath('..\..\..\RadonProject\NewParallelWork\NVIDIA_2DPolarDFT'))
% addpath(genpath('..\..\..\RadonProject\NewParallelWork\NVIDIA_1DFrFT'))
addpath(genpath('..\..\..\RadonProject\NewParallelWork\MEXArrayFireCUDA'))


%% Adding path to Special Correlation On Sphere using SOFT project - My modified form of original and external code
addpath(genpath('..\..\..\MEX_SOFT_Project'))

%% Defining Fixed Rotation Parameters and other initializations
theta_z1 = 5.8*pi/8; % in degrees = 40.5000  >> 0
theta_y = 4.843*pi/8; % in degrees = 63  >> 0
theta_z = 3.84*pi/8; % in degrees = 86.4000  >> 0

% theta_z1 = 3.458*pi/8; % in degrees = 40.5000  >> 0
% theta_y = 5.823*pi/8; % in degrees = 63  >> 0
% theta_z = 2.84*pi/8; % in degrees = 86.4000  >> 0

surfaceVal = 50;
degree = 31;
isPlotting = 1;
scale = .8;
origin=[0 0 0];
vxsize =[1 1 1];

%% Load the original volume data
InDataDir = folderName;

% Without noise
DataFolder = 'benchmark\3DObject_Dataset\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];

inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};

F = length(inFiles);            % Number of Files
Objects_RMSEs       = zeros(1,F);
all_alphas_error    = zeros(1,F);
all_betas_error     = zeros(1,F);
all_gammas_error    = zeros(1,F);
Object_volume_Sizes = zeros(1,F);

viewangle1 = 220;
viewangle2 = 57;

for ii=1:F   % This loop goes through entire Princeton Database, and computes a lot, it takes a lot of time dont run this unless you really want to
%     ii = 74;  %% In the last run these ids had the largest errors , investigated
%     ii = 79;
%     ii = 104;
%     ii = 111;
    %% Loading the Princeton Database Object
    fileName = [InDataDir   DataFolder '\' inFiles(ii).name];
    load(fileName);
    if exist('verticesOrg','var') == 1
        vertices = verticesOrg;
    end
    if exist('facesOrg','var') == 1
        faces = facesOrg;
    end
    [pa,name,ex]=fileparts(fileName);
    
    if(isPlotting)
        %% Before voxelization
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
    end
    %% After Voxelization
    [bim] = verticestovolumefunc(vertices,faces);
    if(isPlotting)
        [vertices, faces] = gen_surf_data(bim,origin,vxsize);
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
    end
    volumeOriginal = surfaceVal *double(bim);
    
    %% The transformation
    t = eye(4);
    t(1:3,1:3)= GetFullRotationMatrixZYZ( theta_z1, theta_y, theta_z );
    tform = affine3d(t);
    
    %% Apply the transformation and display
    volumeRotated = imwarp(volumeOriginal,tform);     % Rotation with pure signal
    [vertices, faces] =  gen_surf_data(logical(volumeRotated/surfaceVal),origin,vxsize);
    
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
    [ final_alpha, final_beta, final_gamma, final_RMSE ] = FindLowestRMSE_Estimate( volumeOriginal, volumeRotated, OutputMatrix, isPlotting );
    
    %% Show Final Transformation 
    
    t = eye(4);
    t(1:3,1:3)= GetFullRotationMatrixZYZ( final_alpha, final_beta, final_gamma );
    tform = affine3d(t);
    volumeRotated_GivenTransformation = imwarp(volumeOriginal,tform);     % Rotation with pure signal
    if(isPlotting)
        %% Plot and Visualize
        [vertices, faces] =  gen_surf_data(logical(volumeOriginal/surfaceVal),origin,vxsize);
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
        
        [vertices, faces] =  gen_surf_data(logical(volumeRotated_GivenTransformation/surfaceVal),origin,vxsize);
        PlotSurface1(vertices,faces);
        view(viewangle1,viewangle2);
    end 
    
    Objects_RMSEs(ii) = final_RMSE;
    all_alphas_error(ii)   = final_alpha - theta_z1;
    all_betas_error(ii)    = final_beta - theta_y;
    all_gammas_error(ii)   = final_gamma - theta_z;
    Object_volume_Sizes(ii) = N;
    
    disp(['Finished Processing, File #', num2str(ii)]);
end

%% Plotting Stuff
figure, plot(1:F, Objects_RMSEs, 'LineWidth', 1)
hold on
plot(1:F, abs(all_alphas_error),'r')
plot(1:F, abs(all_betas_error),'b')
plot(1:F, abs(all_gammas_error),'g')
hold off 
xlabel('Object IDs, from Princeton Shape Database')
ylabel('Error Estimates') 
grid on
legend('RMSE after alignment', '|\alpha - \alpha_{est}| error', '|\beta - \beta_{est}| error', '|\gamma - \gamma_{est}| error')


%% Saving all the objects
save DataComputed\Objects_RMSEs Objects_RMSEs
save DataComputed\all_alphas_error all_alphas_error
save DataComputed\all_betas_error all_betas_error
save DataComputed\all_gammas_error all_gammas_error
save DataComputed\Object_volume_Sizes Object_volume_Sizes
disp('mean'); disp(mean(Objects_RMSEs))
disp('variance'); disp(var(Objects_RMSEs))
disp('mean'); disp(mean(all_alphas_error))
disp('variance'); disp(var(all_alphas_error))
disp('mean'); disp(mean(all_betas_error))
disp('variance'); disp(var(all_betas_error))
disp('mean'); disp(mean(all_gammas_error))
disp('variance'); disp(var(all_gammas_error))
disp('mean'); disp(mean(Object_volume_Sizes))

print -dpdf NoiselessFullComparison