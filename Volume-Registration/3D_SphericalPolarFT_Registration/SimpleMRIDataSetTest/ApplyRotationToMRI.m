close all
clc
clear

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


%% Load the data
% Description :  MR study of head with skull partially removed to reveal brain
% Dimensions  :	 109 slices of 256 x 256 pixels,
%           	 voxel grid is rectangular, and
% 		         X:Y:Z aspect ratio of each voxel is 1:1:2
% Files       :  99 binary files, one file per slice
% File format :  16-bit integers (Mac byte ordering), file contains no header
% Data source :  acquired on a Siemens Magnetom and provided courtesy of
%                Siemens Medical Systems, Inc., Iselin, NJ.  Data edited
%                (skull removed) by Dr. Julian Rosenman, North Carolina
%                Memorial Hospital

% s = load('mri');
% mriVolume = squeeze(s.D);

slicesSize = 99;
N = 69;
mriVolume = zeros(N,N,slicesSize);

mriVolumeNoisy = zeros(N,N,slicesSize);
for s = 1:slicesSize 
    if(s <= 9)
        str = '0';
    else
        str = '';
        
        I = (imresize(flipud(imadjust(imread(strcat('MRIBrainDataSet\mrbrain-16bit0',str,num2str(s),'.tif' )))),[N,N])) ;
        mriVolume(:,:,s) = I;
%         figure, imshow(I); title(num2str(s));
        %% Adding Noise salt pepper
%         J = imnoise(I,'salt & pepper',0.2);
%         mriVolumeNoisy(:,:,s) = J;
%         %% Creating Holes
%         J = I;
%        J(23:39, 43:59) = zeros(length(23:39),length(43:59));
%        mriVolumeNoisy(:,:,s) = J;
% %         figure, imshow(J); title(num2str(s));
        %% PArtial Overlap
% %          %% Creating Holes  Left hole
%         J = I;
%        J((N-1)/2+2:N, 14:N) = zeros(length((N-1)/2+2:N),length(14:N));
%        mriVolumeNoisy(:,:,s) = J;
        %% Creating Holes2 Right Hole
%         J = I;
%         randomInteger = randi(24);
%        J(1:N,(N-1)/2+randomInteger :N ) = zeros(length(1:N),length((N-1)/2+randomInteger :N));
%        mriVolumeNoisy(:,:,s) = J;
       
        %% Creating Holes3 middle
        J = I;
        randomIntegerL = randi(33);
        randomIntegerR = randi(33);
        J((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR,(N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR) = zeros(length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR),length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR));
       mriVolumeNoisy(:,:,s) = J;
%         figure, imshow(J); title(num2str(s));
    end 
end

Visualize3SlicePlanesOfAVolume( mriVolume ); % title('Pure Signal')
  
Visualize3SlicePlanesOfAVolume( mriVolumeNoisy ); % title('Noisy Signal')

%% Volume representation of the data
  
surfaceVal = 15403;
SpecialIsoSurfaceDisplay( mriVolume,surfaceVal );

%% The transformATION , Case 1 , Case 2 , Case 3
theta_z1 = 1.8*pi/8; % in degrees = 40.5000  >> 0
theta_y = 2.8*pi/8; % in degrees = 63  >> 0
theta_z = 3.84*pi/8; % in degrees = 86.4000  >> 0

% theta_x = 2.23*pi/8; % in degrees = 40.5000  >> 0
% theta_y = 1.8*pi/8; % in degrees = 63  >> 0
% theta_z = 3.1*pi/8; % in degrees = 86.4000  >> 0

% theta_x = 4.23*pi/8; % in degrees = 40.5000  >> 0
% theta_y = 3.8*pi/8; % in degrees = 63  >> 0
% theta_z = 3.85*pi/8; % in degrees = 86.4000  >> 0


t = eye(4);
% t(1:3,1:3)= GetFullRotationMatrix( theta_x, theta_y, theta_z );  Not to be used
t(1:3,1:3)= GetFullRotationMatrixZYZ( theta_z1, theta_y, theta_z );
tform = affine3d(t);


%% Apply the transformation
% mriVolumeRotated = imwarp(mriVolume,tform);     % Rotation with pure signal

mriVolumeRotated = imwarp(mriVolumeNoisy,tform);

%% Visualize three slice planes through the center of the transformed volumes.
Visualize3SlicePlanesOfAVolume( mriVolumeRotated ); % title('Rotated Signal')


%% PaddArray To make uniform volumes of N^3 dimensions, N being odd
N = 139;

%% Padding original to match N^3
mriVolume = padarray(mriVolume, [35 35 ],'both');
mriVolume = padarray(mriVolume, [0 0 20],'both');
size(mriVolume)

%% Padding rotated to match N^3  , Case 1 , Case 2 and Case 3
mriVolumeRotated = padarray(mriVolumeRotated, [3 18 3],'both');
mriVolumeRotated = padarray(mriVolumeRotated, [1 0 1],'pre');
size(mriVolumeRotated)
 
% mriVolumeRotated = padarray(mriVolumeRotated, [3 11 0 ],'both');
% mriVolumeRotated = padarray(mriVolumeRotated, [0 1 0],'pre');
% size(mriVolumeRotated)
 
% mriVolumeRotated = padarray(mriVolumeRotated, [15 29 28 ],'both');
% mriVolumeRotated = padarray(mriVolumeRotated, [1 0 0],'pre');
% size(mriVolumeRotated)
 

%% Creating the isotropic low pass filtering in R^3
sigma = .05;
mriVolumeSmooth = imgaussfilt3(mriVolume, sigma);
mriVolumeRotatedSmooth = imgaussfilt3(mriVolumeRotated, sigma);

%% Filterered version plots
Visualize3SlicePlanesOfAVolume( mriVolumeSmooth ); % title('Filtered Padded MRI')

Visualize3SlicePlanesOfAVolume( mriVolumeRotatedSmooth ); % title('Filtered Padded MRI Rotated')

%% Volume Isosurfaces visualization
surfaceVal = 15404;

SpecialIsoSurfaceDisplay( mriVolumeSmooth,surfaceVal );

SpecialIsoSurfaceDisplay( mriVolumeRotatedSmooth,surfaceVal );
view(-185,14);  % Case 1
 
%% Resize for fast computations
ScalingMat = [.45   0   0  0
    0  .45   0  0
    0   0  .45  0
    0   0   0  1];
tformTemp =  affine3d(ScalingMat);
mriVolumeSmooth  = imwarp(mriVolumeSmooth ,tformTemp);
mriVolumeRotatedSmooth = imwarp(mriVolumeRotatedSmooth ,tformTemp);
N = 63;
 
%% Special Function test Final solution

B = N+5; % Must be even
degree = 31;
isPlotting = 1;
OutputMatrix = ComputeSOFTRotation_SphericalPolarFT( mriVolumeSmooth, mriVolumeRotatedSmooth, theta_z1, theta_y, theta_z,  B , degree, isPlotting);

% ComputeRotationSOFTVolumes( mriVolumeSmooth, mriVolumeRotatedSmooth, theta_z1, theta_y, theta_z,  B , degree);

% 
% %% For figures with noisy volume do the following
% 
% print -djpeg 3DNoisyVolumeInputSlices        % because the file size for pdf is too big
% print -dpdf 3DNoisyVolumeSurface
% print -dpdf NoisyVolumeTransReg_1
% print -dpdf NoisyVolumeTransReg_2
% print -dpdf NoisyVolumeTransReg_3
% print -dpdf NoisyVolumeTransReg_4
% print -dpdf NoisyVolumeTransReg_5
% print -dpdf NoisyVolumeTransReg_6
%%
% %% Actual algebraically accurate registration technique, see the paper Volume registration using 3D pseudopolar FFT by Amir Averbuch et.al
% % Computing SphericalGrid for each volume of size zeros(K, M, N+1); Angle phi vs  Polar slices: No of angles theta vs. Radial data
% 
% 
% noOfAnglesTheta               = 2*(N+5);    % Must be even  Special Conversion
% noOfAnglesPhi                 = N+5;      % Must be even
% tic
% mriVolumeSphericalGrid        = abs( VectorizedCompute3DSphericalDFT( (mriVolumeSmooth),  noOfAnglesTheta, noOfAnglesPhi));
% toc
% tic
% mriVolumeRotatedSphericalGrid = abs( VectorizedCompute3DSphericalDFT( (mriVolumeRotatedSmooth),  noOfAnglesTheta, noOfAnglesPhi));
% toc
% 
% %% Special Correlation On Sphere
% addpath(genpath('..\MEX_SOFT_Project'))
% 
% %% Now using the quaternion based representation
% radius = (N-1)/2;  %% Outer most values
% dcVal = mriVolumeSphericalGrid(1,1,(N-1)/2+1);
% degree = 31;  %% Spherical HArmonics coefficient degrees
% LinePoints = -1.4:.02:1.4;
% pointCloudLineX = [LinePoints', zeros(size(LinePoints))',zeros(size(LinePoints))'];
% pointCloudLineY = [zeros(size(LinePoints))',LinePoints',zeros(size(LinePoints))'];
% pointCloudLineZ = [zeros(size(LinePoints))',zeros(size(LinePoints))',LinePoints'];
% pointCloudLine = [pointCloudLineX ;pointCloudLineY; pointCloudLineZ];
% ptCloud = pointCloud(pointCloudLine);
% ptCloudTransformed =  pctransform(ptCloud,tform);
% pointCloudLineTransformed = ptCloudTransformed.Location;
% 
% markerSphereradius = .05;
% lineWidth = 4.5;
% 
% for k =1:radius
%    
%     functionLivingOnSphere1 = GetSphericalImage(  mriVolumeSphericalGrid/dcVal, k );
%      
%     hold on
%     scatter3(pointCloudLine(:,1),pointCloudLine(:,2),pointCloudLine(:,3),'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 0 1])
%     [x,y,z] = sphere(50);
%     x0 = pointCloudLine(end,1);
%     y0 = pointCloudLine(end,2);
%     z0 = pointCloudLine(end,3); 
%     x = x*markerSphereradius + x0;
%     y = y*markerSphereradius + y0;
%     z = z*markerSphereradius + z0;
%     plot3(x,y,z)
%     hold off
%     
%     functionLivingOnSphere2 = GetSphericalImage(  mriVolumeRotatedSphericalGrid/dcVal, k );
%     
%     
%     hold on 
%     scatter3(pointCloudLineTransformed(:,1),pointCloudLineTransformed(:,2),pointCloudLineTransformed(:,3),64,'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 1 0])
%     [x,y,z] = sphere(50);
%     x0 = pointCloudLineTransformed(end,1);
%     y0 = pointCloudLineTransformed(end,2);
%     z0 = pointCloudLineTransformed(end,3);
%     x = x*markerSphereradius + x0;
%     y = y*markerSphereradius + y0;
%     z = z*markerSphereradius + z0;
%     plot3(x,y,z)
%     hold off
%     
% %     previously
% %     I = reshape(functionLivingOnSphere2',[noOfAnglesTheta^2,1]); 
% %     rotatedReal = (reshape(functionLivingOnSphere1',[noOfAnglesTheta^2,1]));
%     
%     I = reshape(functionLivingOnSphere2,[noOfAnglesTheta^2,1]); 
%     rotatedReal = (reshape(functionLivingOnSphere1,[noOfAnglesTheta^2,1]));
%         
% %     [newrotatedReal, newrotatedImag] = MEX_SOFT_Project(1,noOfAnglesTheta/2,noOfAnglesTheta/2,degree,theta_x, theta_y, theta_z,I, zeros(size(I)));   %% Rotating syntax
% %     display_I = imadjust(newrotatedReal);
% %     greyColor = [.7 .7 .7];
% %     [x,y,z]= sphere(30);
% %     figure, 
% %     hold on
% %     warp(x,y,z,(imcomplement(display_I)));
% %     surface(x,y,z,'FaceColor', 'none','EdgeColor',greyColor);
% %     hold off
% %     axis equal;
% %     grid on
% %     xlabel('x')
% %     ylabel('y') 
% %     zlabel('z')
% %     
% %     [alpha, beta, gamma] = MEX_SOFT_Project(0,newrotatedReal, zeros(size(I)),I, zeros(size(I)), noOfAnglesTheta/2,noOfAnglesTheta/2,degree);   %% E syntax
%     
%     [alpha, beta, gamma] = MEX_SOFT_Project(0,rotatedReal, zeros(size(I)),I, zeros(size(I)), noOfAnglesTheta/2,noOfAnglesTheta/2,degree);   %% E syntax
% %     [gamma, beta, alpha] = MEX_SOFT_Project(0,rotatedReal, zeros(size(I)),I, zeros(size(I)), noOfAnglesTheta/2,noOfAnglesTheta/2,degree);   %% E syntax
%     
% %     alpha = 2*pi-alpha;                  %% Resolving ambiguities
% %     beta = pi-beta;
% %     gamma = 2*pi-gamma;
% 
%     R = GetFullRotationMatrixZYZ( alpha, beta, gamma ); 
% %     R = GetFullRotationMatrixZYZ( alpha, beta, gamma ); 
%     
%     t = eye(4);
%     t(1:3,1:3)= R;
%     tform1 = affine3d(t);
%     ptCloudTransformed =  pctransform(ptCloud,tform1);
%     pointCloudLineTransformed = ptCloudTransformed.Location;
%     
%     hold on
%     scatter3(pointCloudLineTransformed(:,1),pointCloudLineTransformed(:,2),pointCloudLineTransformed(:,3),'s','MarkerEdgeColor','k',...
%         'MarkerFaceColor',[1 0 0])
%     [x,y,z] = sphere(50);
%     x0 = pointCloudLineTransformed(end,1);
%     y0 = pointCloudLineTransformed(end,2);
%     z0 = pointCloudLineTransformed(end,3);
%     x = x*markerSphereradius + x0;
%     y = y*markerSphereradius + y0;
%     z = z*markerSphereradius + z0;
%     plot3(x,y,z)
%     hold off
% 
%      
% %     [pointCloud1, sph_vertices1, faces] = GetSphericalMagPointCloud( mriVolumeSphericalGrid/dcVal, k );
% %     %     [pointCloud2, sph_vertices2] = GetSphericalMagPointCloud( mriVolumeRotatedSphericalGrid/dcVal, k );
% %     
% %     pointCloud2 = pctransform(pointCloud(pointCloud1),tform);
% %     [pointCloud2, sph_vertices2, faces ] = SpecialSphericalParametrization( double(pointCloud2.Location), faces );
% %     
% %     SH_Coeffs1 =  GetSphericalHarmonicCoefficients(pointCloud1, sph_vertices1, degree );
% %     SH_Coeffs2 =  GetSphericalHarmonicCoefficients(pointCloud2, sph_vertices2, degree );
% %     
% %     %     [newTransformed,M] = ComputeOptimalRotationQuat(single(pointCloud1),single(pointCloud2));
% %     [newTransformed,M] = ComputeOptimalRotationQuat(SH_Coeffs1,SH_Coeffs2 );              % Comparing SH coefficients for estimating rotations
% %     
% %     %     ptCloudTformed = pctransform(pointCloud(single(pointCloud1)),tform);
% %     %     figure,
%     %     showPointCloud(ptCloudTformed);
%     %     title('Matching Transformed Mag Values');
%     %     xlabel('X')
%     %     ylabel('Y') 
%     %     zlabel('Z')
%     
%     
%     disp(['At radius=', num2str(k)])
%     disp('True transformation')
%     disp(tform.T(1:3,1:3));
%     
%     disp('By SO(3) FFT')
%     disp(R);
%     
%     disp('Acutal')
%     disp([theta_z1,theta_y, theta_z])
%     
%     disp('Estimated')
%     disp([alpha, beta, gamma])
% 
%     disp('Norm absolute error between actual and estimated is = ')
%     disp(norm(abs(tform.T(1:3,1:3))-abs(R)))
%     title(['Registration at spectral radius =',num2str(k)])
%     str = strcat('C:\Users\Alam Abbas\Google Drive\Latex\ImageRegistrationPolarAndSphericalFT\SphericalImageWrapped',num2str(k),'Rotated.pdf');
%     print(gcf,'-dpdf',str);  
% %     disp('By Special quaternion technique')
% %     disp(M'); 
% %     
% %     [theta] = factor_rot_xyz(flipud(R));
% %     disp('Actual')
% %     disp([theta_x, theta_y, theta_z]'*180/pi);
% %     disp('Estimated')
% %     disp(theta);
%     
%     %     if(norm(M'-tform.T) < 10^-2)
%     %         break
%     %     end
% end   
%  
% %% Splitting the whole spherical grid lines from -(N-1)/2:(N-1)/2  to 0:(N-1)/2
% % radialResolution                  = (N-1)/2+1;
% % mriVolumeFullSphericalGrid        = zeros(2*noOfAnglesTheta,2*noOfAnglesPhi,(N-1)/2+1);         % doubled the size with angles from 0 to 360 degrees
% % mriVolumeRotatedFullSphericalGrid = zeros(2*noOfAnglesTheta,2*noOfAnglesPhi,(N-1)/2+1);         % doubled the size with angles from 0 to 360 degrees
% %
% % for m = 1: noOfAnglesTheta
% %     for n = 1: noOfAnglesPhi
% %         %% first volume
% %         radialline = squeeze(mriVolumeSphericalGrid(m,n,:));
% %         firstHalf = radialline(1:(N-1)/2+1);        % Include zero
% %         secondHalf = radialline((N-1)/2+1:N);
% %
% %         mriVolumeFullSphericalGrid (m,n,:) = firstHalf ;
% %         mriVolumeFullSphericalGrid (m+noOfAnglesTheta,n+noOfAnglesPhi,:) = secondHalf;
% %         %% second volume
% %         radialline = squeeze(mriVolumeRotatedSphericalGrid(m,n,:));
% %         firstHalf = radialline(1:(N-1)/2+1);        % Include zero
% %         secondHalf = radialline((N-1)/2+1:N);
% %
% %         mriVolumeRotatedFullSphericalGrid (m,n,:) = firstHalf ;
% %         mriVolumeRotatedFullSphericalGrid (m+noOfAnglesTheta,n+noOfAnglesPhi,:) = secondHalf;
% %     end
% % end
% 
% % xmriVolumeIntegralSphericalGrid = zeros(noOfAnglesTheta,noOfAnglesPhi);
% % mriVolumeRotatedIntegralSphericalGrid = zeros(noOfAnglesTheta,noOfAnglesPhi);
% % absDiffMat = zeros(noOfAnglesTheta,noOfAnglesPhi);
% % radial_values = -(N-1)/2:0; %(N-1)/2;
% % 
% % for m = 1: noOfAnglesTheta
% %     for n = 1: noOfAnglesPhi
% %         line1 = mriVolumeSphericalGrid(m,n,:);
% %         mriVolumeIntegralSphericalGrid(m,n)= trapz(radial_values, line1(1:(N-1)/2+1));
% %         line2 = mriVolumeRotatedSphericalGrid(m,n,:);
% %         mriVolumeRotatedIntegralSphericalGrid(m,n)=  trapz(radial_values, line2(1:(N-1)/2+1));
% %         absDiffMat(m,n) = sum(trapz(line1-line2));
% %         if(m == 34 && n == 56)
% %             figure, plot (squeeze(line1))
% %         end
% %     end
% % end
% % normCorrelation = normxcorr2(mriVolumeIntegralSphericalGrid,mriVolumeRotatedIntegralSphericalGrid);
% % 
% % figure, surf(mriVolumeIntegralSphericalGrid), shading flat
% % view(0,90)
% % xlabel('\theta')
% % ylabel('\phi')
% % 
% % figure, surf(mriVolumeRotatedIntegralSphericalGrid), shading flat
% % view(0,90)
% % xlabel('\theta')
% % ylabel('\phi')
% % 
% % 
% % 
% % [imgIndRows,imgIndCols] = size(mriVolumeRotatedIntegralSphericalGrid);
% % [X,Y,Z] = sphere(imgIndRows);
% % figure,
% % % surface(X,Y,Z,flipud(mriVolumeRotatedIntegralSphericalGrid),...
% % %     'FaceColor','texturemap',...
% % %     'EdgeColor','none',...
% % %     'CDataMapping','direct')
% % % colormap(map)
% % surf(X,Y,Z,flipud(mriVolumeRotatedIntegralSphericalGrid),...
% %     'FaceColor','texturemap',...
% %     'EdgeColor','none',...
% %     'CDataMapping','direct')
% % axis equal
% % view(144,21)
% % grid on
% % rotate3d
% % xlabel('x')
% % ylabel('y')
% % zlabel('z')
% % 
% % 
% % 
% % figure, surf(absDiffMat), shading flat
% % view(0,90)
% % xlabel('\theta')
% % ylabel('\phi')
% % %
% % %
% % % figure, surf(-log(normCorrelation + 20)), shading flat
% % % xlabel('\theta')
% % % ylabel('\phi')
% % 
% % dTheta = pi/noOfAnglesTheta;
% % dPhi  = pi/noOfAnglesPhi;
% % mriVolumeIntegralSphericalGrid = zeros(noOfAnglesTheta,noOfAnglesPhi);
% % mriVolumeRotatedIntegralSphericalGrid = zeros(noOfAnglesTheta,noOfAnglesPhi);
% % 
% % 
% % R = 55; % arbitrary radius
% % 
% % for m = 1: noOfAnglesTheta
% %     for n = 1: noOfAnglesPhi
% %         mriVolumeIntegralSphericalGrid(m,n)= sum(mriVolumeSphericalGrid(m,n,[1:(N-1)/2+1] ));      % Radius 5 information
% %         mriVolumeRotatedIntegralSphericalGrid(m,n)= sum(mriVolumeRotatedSphericalGrid(m,n,[1:(N-1)/2+1] ));
% %     end
% % end
% % 
% % 
% % [theta,phi]= ndgrid(0:dTheta:pi,0:dPhi:pi);
% % 
% % figure,
% % [X,Y,Z] = sph2cart(theta,phi,R );
% % surf(X,Y,Z,mriVolumeIntegralSphericalGrid)
% % xlabel('\theta')
% % ylabel('\phi')
% % axis equal
% % title('Plotting on a sphere: Spectral image 1')
% % colorbar;  daspect([1 1 1]);
% % 
% % figure,
% % [X,Y,Z] = sph2cart(theta,phi,R );
% % surf(X,Y,Z,mriVolumeRotatedIntegralSphericalGrid)
% % xlabel('\theta')
% % ylabel('\phi')
% % axis equal
% % title('Plotting on a sphere: Spectral image 2')
% % colorbar;  daspect([1 1 1]);
% % 
% % normCorrelation = normxcorr2(mriVolumeIntegralSphericalGrid,mriVolumeRotatedIntegralSphericalGrid);
% % figure, surf(normCorrelation), shading flat
% % view(0,90)
% % xlabel('\theta')
% % ylabel('\phi')