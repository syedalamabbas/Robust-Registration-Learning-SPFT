%% Serialization test scipt
% Written on 9/18/2022 by Alam Abbas Syed

clear
clc
close all


%% Fetch a set of spherical images from SPFT
priorfolder1 = 'extracted_multispectral_sp_images\';
priorfolder2 = 'extracted_multispectral_sp_noise_images\';
priorfolder3 = 'extracted_rscnn_sp_images\';
priorfolder4 = 'extracted_rscnn_sp_noise_images\';
f_name = '\airplane\test\airplane_0627_r';


for k =1:4
    str_k = num2str(k);
    sp_img_1d = csvread(strcat(strcat(priorfolder1 , f_name), str_k, '.csv'));
    figure,
    imagesc(reshape(sp_img_1d, [64 64]));
%     colorbar;
    axis off
    sp_img_1d = csvread(strcat(strcat(priorfolder2 , f_name), str_k, '.csv'));
    figure,
    imagesc(reshape(sp_img_1d, [64 64]));
%     colorbar;
    axis off
    sp_img_1d = csvread(strcat(strcat(priorfolder3 , f_name), str_k, '.csv'));
    figure,
    imagesc(reshape(sp_img_1d, [64 64]));
%     colorbar;
    axis off
    sp_img_1d = csvread(strcat(strcat(priorfolder4 , f_name), str_k, '.csv'));
    figure,
    imagesc(reshape(sp_img_1d, [64 64]));
%     colorbar;
    axis off
end

%%
tic
N = 65; % Output is N x N x N volume
voxels = randn(N, N, N);

% Serialize
strFileToSave = 'SerializeDeserialize.csv';
csvwrite(strFileToSave ,voxels(:));

% Deserialize
deserialized_voxels_1D = csvread(strFileToSave);
deserialized_voxels = reshape(deserialized_voxels_1D, [N, N, N]);

% Compare
disp(['Sum difference between serialized and deserialized voxels is ', num2str(sum(sum(sum(voxels - deserialized_voxels))))])
% Should be exactly or close to zero

delete(strFileToSave)
toc
