function [ OutputMatrix ] = ComputeSOFTRotationPairs( VolumeSphericalGrid, VolumeRotatedSphericalGrid, theta_z1, theta_y, theta_z, B , degree, isPlotting )
%COMPUTESOFTROTATIONPAIRS Computes the rotation estimation using Spherical
%Harmonics
% VolumeSphericalGrid        - The M x K x N volume that is the reference 
% VolumeRotatedSphericalGrid - The second M x K x N volume that is the rotated volume w.r.t first volume
% theta_z1, theta_y, theta_z - These angles are the true transformations only needed for plotting purposes
% B                          - Bandwidth for the SOFT technique
% degree                     - Spherical Harmonics coefficient degrees for SOFT technique
% isPlotting                 - Plotting the solution = 1 or not = 0
% OutputMatrix               - Spectral Radius vs [alpha, beta, gamma]
%====================================================================
% Written on June 22nd, 2016 by Syed Alam Abbas.
% Revised on March 12th, 2022 by Syed Alam Abbas.
%====================================================================

%% The two volumes must be of the same size, the bandwidth B must be even Check inputs
noOfAnglesTheta               = 2*B;    % Must be even  Special Conversion and the ratio
noOfAnglesPhi                 = B;      % Must be even
[~,~,N] = size(VolumeSphericalGrid);
%% Creating Transform matrix
t = eye(4);
t(1:3,1:3)= GetFullRotationMatrixZYZ( theta_z1, theta_y, theta_z );
tform = affine3d(t);

radius = 11; (N-1)/2;  %% Outer most values
dcVal = VolumeSphericalGrid(1,1,(N-1)/2+1);

OutputMatrix = zeros(radius,3);    % Spectral Radius x [alpha, beta, gamma]

if(isPlotting) 
    LinePoints = -1.4:.02:1.4;                % Points around and covering the unit sphere
    pointCloudLineX = [LinePoints', zeros(size(LinePoints))',zeros(size(LinePoints))'];
    pointCloudLineY = [zeros(size(LinePoints))',LinePoints',zeros(size(LinePoints))'];
    pointCloudLineZ = [zeros(size(LinePoints))',zeros(size(LinePoints))',LinePoints'];
    pointCloudLine = [pointCloudLineX ;pointCloudLineY; pointCloudLineZ];
    ptCloud = pointCloud(pointCloudLine);
    ptCloudTransformed =  pctransform(ptCloud,tform);
    pointCloudLineTransformed = ptCloudTransformed.Location;
end

for k =1:radius
    
    %% Get Spherical Images at the spectral radius = k
    functionLivingOnSphere1 = GetSphericalImage(  VolumeSphericalGrid/dcVal, k , isPlotting);
    if(isPlotting)
       PlotMarkerLinesAtOrigin(pointCloudLine, 36, [0 0 1], 'o');   % Plotting marker lines over the images obtained from the above function 
    end
    functionLivingOnSphere2 = GetSphericalImage(  VolumeRotatedSphericalGrid/dcVal, k , isPlotting);
    if(isPlotting)
       PlotMarkerLinesAtOrigin(pointCloudLineTransformed,64,[0 1 0],'o');  % Plotting marker lines over the images obtained from the above function 
    end
    
    %%  Reshape Images into  noOfAnglesTheta x noOfAnglesTheta matrix
    I = reshape(functionLivingOnSphere2',[noOfAnglesTheta^2,1]);                % Not sure why do we need a transpose
    complexPart_SOFT =  zeros(size(I));
    rotatedReal = (reshape(functionLivingOnSphere1',[noOfAnglesTheta^2,1]));    % Not sure why do we need a transpose
    
    %% Execute SOFT Technique on the two obtained spherical images
    [alpha, beta, gamma] = MEX_SOFT_Project(0,rotatedReal, complexPart_SOFT ,I, complexPart_SOFT, B,B,degree);   %% E syntax
    R = GetFullRotationMatrixZYZ( alpha, beta, gamma );
    if(isPlotting)
        t = eye(4);
        t(1:3,1:3)= R;
        tformRotated = affine3d(t);
        ptCloudTransformed2 =  pctransform(ptCloud,tformRotated);
        pointCloudLineTransformed2 = ptCloudTransformed2.Location;
        PlotMarkerLinesAtOrigin(pointCloudLineTransformed2,36,[1 0 0],'s');     % Plotting marker lines over the images obtained from the above function 
        title(['Registration at spectral radius =',num2str(k)])
    end
    
    %% Displaying some information about the computed solution
    disp(['At radius=', num2str(k)])
    disp('True transformation')
    disp(tform.T(1:3,1:3));
    
    disp('By my SO(3) FFT')
    disp(R);
    
    disp('Acutal')
    disp([theta_z1,theta_y, theta_z])
      
    disp('Estimated by my SO(3) FFT')
    disp([alpha, beta, gamma])
    
    disp('Norm absolute error between actual and estimated is = ')
    disp(norm(abs(tform.T(1:3,1:3))-abs(R)))
    
    %% Publish the results possibly in a file or give an output
    OutputMatrix(k,:) = [alpha,beta,gamma];
end 

    function PlotMarkerLinesAtOrigin(pointCloudLocation, size, faceColor, shape)
        hold on 
        scatter3(pointCloudLocation(:,1),pointCloudLocation(:,2),pointCloudLocation(:,3),size,shape, 'MarkerEdgeColor','k',...
            'MarkerFaceColor',faceColor)
        [x,y,z] = sphere(50);               % Marker sphere at the end of the point cloud sequence
        markerSphereradius = .05;
        x0 = pointCloudLocation(end,1);
        y0 = pointCloudLocation(end,2);
        z0 = pointCloudLocation(end,3);
        x = x*markerSphereradius + x0;
        y = y*markerSphereradius + y0;
        z = z*markerSphereradius + z0;
        plot3(x,y,z)
        hold off
    end
end

