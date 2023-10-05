%% Preprocessing scipt
% Written on 9/17/2022 by Alam Abbas Syed

clear
clc
close all

%% Initialization
addpath(genpath('archive'))
addpath(genpath('../Common'))
metadata_modelnet40 = readtable('metadata_modelnet40.csv');
N = 65; % Output is N x N x N volume
sz = [N N N];

%% Read files
preFolder = 'archive/ModelNet40/';
[table_rows, table_cols] = size(metadata_modelnet40);
% Chair = 3931
% Airplane = 1
plot_ids = [ randi(table_rows), randi(table_rows), randi(table_rows), randi(table_rows), randi(table_rows) ];

%% Display random objects
for k=1:length(plot_ids)
    strFilePath = string(metadata_modelnet40{plot_ids(k), 4});
    disp(['Object id=', num2str(plot_ids(k))])
    disp([ ' path=', strFilePath ]);
    
    [vertices,faces] = read_off( preFolder + strFilePath );
    fvc3.vertices=vertices';
    fvc3.faces=faces';
    
    % Normalizing input vertices
    [r, c] = size(fvc3.vertices);
    avgx = sum(fvc3.vertices(:, 1))/r;
    avgy = sum(fvc3.vertices(:, 2))/r;
    avgz = sum(fvc3.vertices(:, 3))/r;
    fvc3.vertices(:, 1) = fvc3.vertices(:, 1) - avgx;
    fvc3.vertices(:, 2) = fvc3.vertices(:, 2) - avgy;
    fvc3.vertices(:, 3) = fvc3.vertices(:, 3) - avgz;
    
    % Before Voxelization
    PlotSurface1(fvc3.vertices,fvc3.faces);
    
    % Voxelize
    voxels_volume=polygon2voxel(fvc3, sz,'auto'); % Logical volume
    
    % After Voxelization
    origin =[0 0 0];
    vxsize =[1 1 1];
    [vox_vertices, vox_faces] =  gen_surf_data(voxels_volume,origin,vxsize);
    PlotSurface1(vox_vertices, vox_faces);
    
    % Adding rotations, missing, outliers, noise for ablation study
%     z_rot = rand() * 2 * pi; % Rotate vertices along Z axis (or Azimuth) between 0 to 2*pi
%     Rz = [ cos(z_rot) -sin(z_rot)       0; ...
%        sin(z_rot)  cos(z_rot)       0; ...
%             0      0        1];
%     vox_vertices = vox_vertices * Rz; %Clockwise rotation exists
%     out_vox_vertices = round(AddAWGN_XYZPoints( vox_vertices,  40 )); % SNRdB = 10,20,30,40
    
%     out_vox_vertices=round(outliers(vox_vertices',.5,[1 N]))';
%     out_vox_vertices = round(AddAWGN_XYZPoints( vox_vertices,  20 )); % SNRdB = 10,20,30,40
    out_vox_vertices=(missing_points(vox_vertices',.9))';
      
    mod_vol = logical(zeros(sz));
    I1 = out_vox_vertices(:,1);
    I1(I1 <= 0) = 1;
    I1(I1 >= 65) = 65;
    I2 = out_vox_vertices(:,2);
    I2(I2 <= 0) = 1;
    I2(I2 >= 65) = 65;
    I3 = out_vox_vertices(:,3);
    I3(I3 <= 0) = 1;
    I3(I3 >= 65) = 65;
    ind = sub2ind(sz,I1,I2,I3);
    mod_vol(ind) = 1;
    
    figure, volshow(voxels_volume)
    figure, pcshow(pointCloud(vox_vertices))
    axis tight
    figure, pcshow(pointCloud(out_vox_vertices))
    axis tight
    figure, volshow(mod_vol)
end

%% Process entire dataset and save processed volumes
close all
tic
disp('Completed processing ModelNet40 dataset in ...')
% processedFolder = 'extracted_volumes\'; % Clean
% processedFolder = 'extracted_volumes_outliers\';
% processedFolder = 'extracted_volumes_noise\';
% processedFolder = 'extracted_volumes_missing\';
processedFolder = 'extracted_volumes_rot_noise\';
if ~exist(processedFolder, 'dir')
    mkdir(processedFolder)
end

fvc3_array = cell(table_rows, 1);
parfor row=1:table_rows
    % Read off file
    strFilePath = string(metadata_modelnet40{row, 4});
    [vertices,faces] = read_off( preFolder + strFilePath);
    fvc3_array{row}.vertices=vertices';
    fvc3_array{row}.faces=faces';
    
    % Normalize input vertices
    [r, c] = size(fvc3_array{row}.vertices);
    avgx = sum(fvc3_array{row}.vertices(:, 1))/r;
    avgy = sum(fvc3_array{row}.vertices(:, 2))/r;
    avgz = sum(fvc3_array{row}.vertices(:, 3))/r;
    fvc3_array{row}.vertices(:, 1) = fvc3_array{row}.vertices(:, 1) - avgx;
    fvc3_array{row}.vertices(:, 2) = fvc3_array{row}.vertices(:, 2) - avgy;
    fvc3_array{row}.vertices(:, 3) = fvc3_array{row}.vertices(:, 3) - avgz;
    
    % Voxelize
    voxels_volume = polygon2voxel(fvc3_array{row},[N  N  N],'auto');
    
    % For ablation study
    [vox_vertices, vox_faces] =  gen_surf_data(voxels_volume,[0 0 0],[1 1 1]);
    
    z_rot = rand() * 2 * pi; % Rotate vertices along Z axis (or Azimuth) between 0 to 2*pi
    Rz = [ cos(z_rot) -sin(z_rot)       0; ...
       sin(z_rot)  cos(z_rot)       0; ...
            0      0        1];
    vox_vertices = vox_vertices * Rz; %Clockwise rotation exists
    out_vox_vertices = round(AddAWGN_XYZPoints( vox_vertices,  40 )); % SNRdB = 10,20,30,40
    
%     out_vox_vertices=round(outliers(vox_vertices',.5,[1 N]))'; % 50% outliers
%     out_vox_vertices = round(AddAWGN_XYZPoints( vox_vertices,  20 )); % SNRdB = 10,20,30,40
%     out_vox_vertices=(missing_points(vox_vertices',.9))';
      
    mod_vol = logical(zeros(sz));
    I1 = out_vox_vertices(:,1);
    I1(I1 <= 0) = 1;
    I1(I1 >= N) = N;
    I2 = out_vox_vertices(:,2);
    I2(I2 <= 0) = 1;
    I2(I2 >= N) = N;
    I3 = out_vox_vertices(:,3);
    I3(I3 <= 0) = 1;
    I3(I3 >= N) = N;
    ind = sub2ind(sz,I1,I2,I3);
    mod_vol(ind) = 1;

    % Modifying the volume
    voxels_volume = mod_vol;
    
    % Save volume
    strObj = string(metadata_modelnet40{row,1});
    strClass = string(metadata_modelnet40{row,2});
    strSplit = string(metadata_modelnet40{row,3});
    
    if ~exist(strcat(processedFolder, strClass), 'dir')
        mkdir(strcat(processedFolder, strClass))
    end
    if ~exist(strcat(processedFolder, strClass, '\', strSplit), 'dir')
        mkdir(strcat(processedFolder, strClass, '\', strSplit))
    end
    
    strFileToSave = strcat(processedFolder, strClass, '\', strSplit, '\', strObj, '.csv' );
    csvwrite(strFileToSave, voxels_volume(:));
end
toc
% Elapsed time is 15086.277750 seconds. (That is roughly 4.19 hrs on my machine!)
% Elapsed time is 34511.452947 seconds. 9.58 hrs on my laptop, file saves are terrible
% Outliers: Elapsed time is 31462.507162 seconds.
% Noise: Elapsed time is 29880.694776 seconds.
% Missing: Elapsed time is 29607.915074 seconds.
% Rot + Noise: Elapsed time is 31235.265323 seconds.