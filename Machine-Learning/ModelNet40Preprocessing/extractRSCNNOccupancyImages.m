%% Preprocessing scipt
% Written on 9/24/2022 by Alam Abbas Syed
% Edited on 9/27/2022

clear
clc
close all

%% Initialization
addpath(genpath('archive'))
addpath(genpath('../Common'))
metadata_modelnet40 = readtable('metadata_modelnet40.csv');
N = 65; % Output is N x N x N volume
nR = 7; % Number of radii for spatial occupancy spherical grid
%% Read files
preFolder = 'archive/ModelNet40/';
[table_rows, table_cols] = size(metadata_modelnet40);

%% Process entire dataset and save processed spatial occupancy grid images
tic
disp('Completed processing ModelNet40 dataset for R-SCNN images ...')
% rscnn_sp_imgs_Folder = 'extracted_rscnn_sp_images\';
% rscnn_sp_imgs_Folder = 'extracted_rscnn_sp_noise_images\';
rscnn_sp_imgs_Folder = 'extracted_rscnn_sp_rot_noise_images\';
if ~exist(rscnn_sp_imgs_Folder, 'dir')
    mkdir(rscnn_sp_imgs_Folder)
end

fvc3_array = cell(table_rows, 1);
for row=1:table_rows
    strObj = string(metadata_modelnet40{row,1});
    strClass = string(metadata_modelnet40{row,2});
    strSplit = string(metadata_modelnet40{row,3});
    
    
    % Read off file
    strFilePath = string(metadata_modelnet40{row, 4});
    %     disp(['At row = ', num2str(row) , ' path=' , strFilePath])
    [vertices, faces] = read_off( preFolder + strFilePath);
    
    fvc3_array{row}.vertices=vertices';
    fvc3_array{row}.faces=faces';
    
    % Normalize input vertices
    [r, ~] = size(fvc3_array{row}.vertices);
    avgx = sum(fvc3_array{row}.vertices(:, 1))/r;
    avgy = sum(fvc3_array{row}.vertices(:, 2))/r;
    avgz = sum(fvc3_array{row}.vertices(:, 3))/r;
    fvc3_array{row}.vertices(:, 1) = fvc3_array{row}.vertices(:, 1) - avgx;
    fvc3_array{row}.vertices(:, 2) = fvc3_array{row}.vertices(:, 2) - avgy;
    fvc3_array{row}.vertices(:, 3) = fvc3_array{row}.vertices(:, 3) - avgz;
    
    % Voxelize
    voxels_volume = polygon2voxel(fvc3_array{row},[N  N  N],'auto');
    [vertices, faces] = gen_surf_data(voxels_volume,[0 0 0],N);
    PlotSurface1(vertices, faces);
    
    % Adding noise and all
    vox_vertices = vertices;
%     z_rot = rand() * 2 * pi; % Rotate vertices along Z axis (or Azimuth) between 0 to 2*pi
%     Rz = [ cos(z_rot) -sin(z_rot)       0; ...
%        sin(z_rot)  cos(z_rot)       0; ...
%             0      0        1];
%     vox_vertices = vox_vertices * Rz; %Clockwise rotation exists
%     out_vox_vertices = round(AddAWGN_XYZPoints( vox_vertices,  40 ));
    
%      out_vox_vertices=round(outliers(vox_vertices',.5,[1 N]))'; % 50% outliers
    out_vox_vertices = round(AddAWGN_XYZPoints( vox_vertices,  20 )); % SNRdB = 10,20,30,40
%     out_vox_vertices=(missing_points(vox_vertices',.9))';

%     out_vox_vertices = round(AddAWGN_XYZPoints( vertices,  20 ));
    vertices = out_vox_vertices;
    PlotSurface1(vertices, faces);
    
    vertices = vertices';
    
    % Normalize scale to fit points in a unit sphere
    [~, c] = size(vertices);
    centroid_c = [sum(vertices(1, :)) sum(vertices(2, :)) sum(vertices(3, :))]/c;
    vertices(1, :) = vertices(1,:) - centroid_c(1);
    vertices(2, :) = vertices(2,:) - centroid_c(2);
    vertices(3, :) = vertices(3,:) - centroid_c(3);
    furthest_distance = max(sqrt(sum(abs(vertices(3, :)).*2)));
    vertices = vertices / furthest_distance;
    
    points = vertices;
    ptCloud=pointCloud(points');
    pcshow(ptCloud);
    normals_c = pcnormals(ptCloud,9)';
    
    v = points./sum(points.^2,1).^.5;     % find point-000 vector
    thetas = acos(abs(sum(normals_c.*v)));   % angle ptw two vectors
    
    [theta, rho, z]= cart2sph(points(1,:),points(2,:),points(3,:));
    points = [theta;rho;z];
    [rr, thetas] = sphvoxels(points,thetas,(N-1),nR);
    
    r = permute(rr, [3, 1,2]);
    thetas = permute(thetas, [3,1,2]);
    
    if(strcmp(strClass,'airplane'))
        pcshow(ptCloud)
        for n=1:nR
            figure, imagesc(squeeze(r(n, :, :)))
            colorbar
        end
    end
    
    % Save images
    if ~exist(strcat(rscnn_sp_imgs_Folder, strClass), 'dir')
        mkdir(strcat(rscnn_sp_imgs_Folder, strClass))
    end
    if ~exist(strcat(rscnn_sp_imgs_Folder, strClass, '\', strSplit), 'dir')
        mkdir(strcat(rscnn_sp_imgs_Folder, strClass, '\', strSplit))
    end
    
    for n =1:nR
        ImageOnSphere = squeeze(r(n, :, :));
        strFileToSave = strcat(rscnn_sp_imgs_Folder, strClass, '\', strSplit, '\', strObj, '_r', num2str(n), '.csv' );
        csvwrite(strFileToSave, ImageOnSphere(:));
    end
end
toc
% Need a better laptop for this...
% Elapsed time is 8514.404728 seconds.
% Noise: Elapsed time is 7477.445813 seconds.