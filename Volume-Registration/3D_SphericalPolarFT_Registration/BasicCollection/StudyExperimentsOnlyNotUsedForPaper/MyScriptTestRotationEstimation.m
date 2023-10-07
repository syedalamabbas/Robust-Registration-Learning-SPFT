
close all
clc

%% Adding path to the root folder of this project i.e. 3D_SphericalPolarFT_Registration
addpath(genpath('..\'))

%% Original
ptCloud = pcread('teapot.ply');
figure
showPointCloud(ptCloud);
title('Teapot');
xlabel('X')
ylabel('Y')
zlabel('Z')


%% Transformation
delay = .05;
figure,
alpha = 3.875*pi/6;
beta = 2.965*pi/6;
gamma = 3.6558*pi/6;
R = GetFullRotationMatrix( alpha, beta, gamma );
t_Vector = [0  0 0];
A = eye(4);
A(1:3,1:3)= R;
A(4,1:3)= t_Vector;

tform1 = affine3d(A);
ptCloudTformed = pctransform(ptCloud,tform1);
showPointCloud(ptCloudTformed);
title('Transformed Teapot');
xlabel('X') 
ylabel('Y') 
zlabel('Z')

%% Special quaternion based optimal rotation angle estimation

[newTransformed,M] = ComputeOptimalRotationQuat(ptCloud.Location,ptCloudTformed.Location);
% degree= 15;
% [ sph_vertices1 ] = SpecialSphericalParametrization(double( ptCloud.Location) );
% [ sph_vertices2 ] = SpecialSphericalParametrization( double( ptCloudTformed.Location) );
% 
% 
% SH_Coeffs1 =  GetSphericalHarmonicCoefficients(pointCloud1, sph_vertices1, degree );
% SH_Coeffs2 =  GetSphericalHarmonicCoefficients(pointCloud2, sph_vertices2, degree );
% 
% %     [newTransformed,M] = ComputeOptimalRotationQuat(single(pointCloud1),single(pointCloud2));
% [newTransformed,M] = ComputeOptimalRotationQuat(SH_Coeffs1,SH_Coeffs2 );              % Comparing SH coefficients for estimating rotations


%% Registration by ICP
tform = pcregrigid(ptCloudTformed,ptCloud,'Extrapolate',true);

%% Final result comparison
disp('True transformation') 
disp(tform1.T);

disp('By Special quaternion technique')
disp(M');

tform2 = invert(tform);
disp('Transformation by ICP algorithm')
disp(tform2.T);
