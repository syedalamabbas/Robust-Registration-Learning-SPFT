clear 
clc
close all

%% Initialize
addpath(genpath('archive'))
addpath(genpath('../Common'))
% addpath(genpath('extraced_spherical_grids'))
N = 65;
B = (N-1)/2;
noOfAnglesTheta = 2*B;    % Must be even  Special Conversion and the ratio
noOfAnglesPhi   = B;      % Must be even

R = 11; % Number of radii for extracting spectral signatures
%% Load and Display spectral signatures for airplane_0627 model
strFileToLoadVoxels = 'airplane_0627.csv';
voxels_volume_1D = csvread(strFileToLoadVoxels);
voxels_volume = reshape(voxels_volume_1D, [N, N, N]);

% Visualize Voxelization as Surface/Mesh
origin =[0 0 0];
vxsize =[1 1 1];
[vox_vertices, vox_faces] =  gen_surf_data(voxels_volume,origin,vxsize);
PlotSurface1(vox_vertices, vox_faces);

% Visualize Voxelization as Point Cloud
figure, 
pcshow(vox_vertices)
xlabel('x')
ylabel('y')
zlabel('z')
% view(45, 0)

% Spherical grid
strFileToLoadSpGrid = 'airplane_0627_sp.csv';
deserialized_spgrid_1D = csvread(strFileToLoadSpGrid);
deserialized_spgrid = reshape(deserialized_spgrid_1D, [noOfAnglesTheta, noOfAnglesPhi, N]);
dcVal = deserialized_spgrid(1,1,(N-1)/2+1);
deserialized_spgrid = deserialized_spgrid / dcVal;

isPlotting = 1;
isRemoveAxesLabels = 1;
for r =1:R
    ImageOnSphere = GetSphericalImage(deserialized_spgrid, r , isPlotting, isRemoveAxesLabels);
end

%% Extracting multispectral spherical images
close all
tic
disp('Completed extraction of multispectral images from ModelNet40 dataset in ...')
metadata_modelnet40 = readtable('metadata_modelnet40.csv');
% spectral_imgs_Folder = 'extracted_multispectral_sp_images\';
% sp_grids_folder = 'extracted_spherical_grids\';
% spectral_imgs_Folder = 'extracted_multispectral_sp_noise_images\';
% sp_grids_folder = 'extracted_spherical_grids_noise\';
% spectral_imgs_Folder = 'extracted_multispectral_sp_missing_images\';
% sp_grids_folder = 'extracted_spherical_grids_missing\'; 
% spectral_imgs_Folder = 'extracted_multispectral_sp_outliers_images\';
% sp_grids_folder = 'extracted_spherical_grids_outliers\'; 
spectral_imgs_Folder = 'extracted_multispectral_sp_rot_noise_images\';
sp_grids_folder = 'extracted_spherical_grids_rot_noise\'; 

if ~exist(spectral_imgs_Folder, 'dir')
    mkdir(spectral_imgs_Folder)
end

[table_rows, table_cols] = size(metadata_modelnet40);
parfor row=1:table_rows
    % Read Meta file info
    strObj = string(metadata_modelnet40{row,1});
    strClass = string(metadata_modelnet40{row,2});
    strSplit = string(metadata_modelnet40{row,3});
    
    % Read spectral full grid
    strFileToLoadSpGrid = strcat(sp_grids_folder, strClass, '\', strSplit, '\', strObj, '_sp', '.csv' );
    deserialized_spgrid_1D = csvread(strFileToLoadSpGrid);
    deserialized_spgrid = reshape(deserialized_spgrid_1D, [noOfAnglesTheta, noOfAnglesPhi, N]);
    dcVal = deserialized_spgrid(1,1,(N-1)/2+1);
    deserialized_spgrid = deserialized_spgrid / dcVal;

    % Save normalized image arrays
    if ~exist(strcat(spectral_imgs_Folder, strClass), 'dir')
        mkdir(strcat(spectral_imgs_Folder, strClass))
    end
    if ~exist(strcat(spectral_imgs_Folder, strClass, '\', strSplit), 'dir')
        mkdir(strcat(spectral_imgs_Folder, strClass, '\', strSplit))
    end
   
    for r =1:R
        ImageOnSphere = GetSphericalImage(deserialized_spgrid, r , 0, 0);
        strFileToSave = strcat(spectral_imgs_Folder, strClass, '\', strSplit, '\', strObj, '_r', num2str(r), '.csv' );
        csvwrite(strFileToSave, ImageOnSphere(:));
    end
end
toc
%Elapsed time is 8259.476413 seconds.