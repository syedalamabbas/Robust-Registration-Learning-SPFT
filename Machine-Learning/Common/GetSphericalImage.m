function [ ImageOnSphere ] = GetSphericalImage( VolumeMagSphericalGrid, radius , isPlotting, isRemoveAxesLabels)
%% Written by Alam Abbas Syed, 2018 - edited on 9/20/2022
%GETSPHERICALMAGPOINTCLOUD computes the point cloud based on intensities on
%the sphere at a specific radius given as an input
% Example call : VolumeMagSphericalGrid = rand(100,100,99);
% Then:          GetSphericalMagPointCloud( VolumeMagSphericalGrid, 45 )
% -- get all points at radius 45

[ noOfAnglesTheta, noOfAnglesPhi, sizeX ] = size(VolumeMagSphericalGrid);
N = sizeX -1;                          % N is even
P = noOfAnglesTheta;                   % P is also even it is for theta
Q = noOfAnglesPhi;                     % Q is also even it is for phi

deltaTheta = 180/P;                    % Angular sampling rate theta
anglesTheta = 0:deltaTheta:180-deltaTheta;

deltaPhi = 180/Q;                      % Angular sampling rate phi
anglesPhi = 0:deltaPhi:180-deltaPhi;

gridSpacing =  -N/2:N/2;
sphericalVertices = zeros(P*Q*2,3);      % Total number of points on a sphere

% Some more special variables
functionLivingOnSphere = zeros(P-1,2*Q);    % P == 2*Q for an equiangular square image on a sphere used for correlation
B = Q;
% Specialized Indexes
ImageOnSphere = zeros(2*B-1,2*B);     % 4*B^2

Quad2IndexesHorizontal = B:-1:1;
Quad2IndexesVertical = B:-1:1;

Quad1IndexesHorizontal = 2*B:-1:B+1;
Quad1IndexesVertical = 1:1:B+1;

Quad4IndexesHorizontal = 2*B:-1:B+1;
Quad4IndexesVertical = B+1:1:2*B;

Quad3IndexesHorizontal = B:-1:1;
Quad3IndexesVertical = 2*B:-1:B+1;

count = 1;
amplitudeOffset = 0;.3;
for p = 1:P
    angleTheta = anglesTheta(p);
    for q = 1:Q
        anglePhi = anglesPhi (q);
        for n = 1:N+1
            %             if(gridSpacing(n) == radius)
            %                 functionLivingOnSphere(p,q) = VolumeMagSphericalGrid(p,q,n)+.3;
            %             end
            %             if(-gridSpacing(n) == radius)
            %                 functionLivingOnSphere(p,q+Q) = VolumeMagSphericalGrid(p,q,n)+.3;
            %             end
            %             if(abs(gridSpacing(n)) == radius )
            %                 %% For plotting the intensity based object
            %                 magValue = sign(gridSpacing(n))* (VolumeMagSphericalGrid(p,q,n)+.3);  %% for actual using the value
            %                 desiredPoint = [magValue*cosd(angleTheta)*cosd(anglePhi),magValue*cosd(angleTheta)*sind(anglePhi),magValue*sind(angleTheta)];
            %                 desiredPointCloud(count,:) = desiredPoint;
            %                 functionValueOnVertex(count) = VolumeMagSphericalGrid(p,q,n)+.3;
            %                 %                 %% For plotting the sphere vertices
            %                 magValue = sign(gridSpacing(n))*radius;
            %                 pointOnSphere = [magValue*cosd(angleTheta)*cosd(anglePhi),magValue*cosd(angleTheta)*sind(anglePhi),magValue*sind(angleTheta)];
            %                 sphericalVertices(count,:) = pointOnSphere;
            %                 %% Increment for looping
            %                 count = count +1;
            %             end
            
            if(+(gridSpacing(n)) == radius )
                magValue = sign(gridSpacing(n))* (VolumeMagSphericalGrid(p,q,n)+amplitudeOffset );
                pointOnSphere = [magValue*cosd(angleTheta)*cosd(anglePhi),magValue*cosd(angleTheta)*sind(anglePhi),magValue*sind(angleTheta)];
                sphericalVertices(count,:) = pointOnSphere;
               
                
                %% Special Function update
                if(p <= Q)
                    ImageOnSphere(Quad2IndexesVertical(p),Quad2IndexesHorizontal(q)) = sign(gridSpacing(n))*magValue;
                end
                if(p > Q+1 && p <= P)
                    ImageOnSphere(Quad1IndexesVertical(p-Q-1),Quad1IndexesHorizontal(q)) =sign(gridSpacing(n))* magValue;
                end
                
                %% Increment for looping
                count = count +1;
            end
            
            if(-(gridSpacing(n)) == radius )
                magValue = sign(gridSpacing(n))* (VolumeMagSphericalGrid(p,q,n)+amplitudeOffset );
                pointOnSphere = [magValue*cosd(angleTheta)*cosd(anglePhi),magValue*cosd(angleTheta)*sind(anglePhi),magValue*sind(angleTheta)];
                sphericalVertices(count,:) = pointOnSphere;
             
                
                %% Special Function update
                if(p <= Q)
                    ImageOnSphere(Quad4IndexesVertical(p)-1,Quad4IndexesHorizontal(q)) = sign(gridSpacing(n))*magValue;
                end
                if(p > Q+1 && p <= P)
                    ImageOnSphere(Quad3IndexesVertical(p-Q-1)-1,Quad3IndexesHorizontal(q)) =sign(gridSpacing(n))* magValue;
                end
                
                %% Increment for looping
                count = count +1;
            end
        end
    end
end

% figure,
% view(-63,18)
% axis equal
% axis ([-1.5 1.5 -1.5 1.5 -1.5 1.5])
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid on
% hold on
% for l = 1: count
%     plot3(sphericalVertices(l,1),sphericalVertices(l,2),sphericalVertices(l,3) ,'.');
%     drawnow
% end
% hold off

% %% wrap it in on sphere
% display_I = imadjust(functionLivingOnSphere);
% greyColor = [.7 .7 .7];
% [x,y,z]= sphere(30);
% figure,
% hold on
% warp(x,y,z,(imcomplement(display_I)));
% surface(x,y,z,'FaceColor', 'none','EdgeColor',greyColor);
% hold off
% axis equal;
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('z')

%% wrap it in on sphere
ImageOnSphere = imresize( ImageOnSphere, [P ,P]);
    
if (isPlotting)
    figure, imagesc(ImageOnSphere);
    if(isRemoveAxesLabels)
        axis off
    else
        colorbar
    end

%     % display_I = imadjust(reshape(functionValueOnVertex,[P, P]));
%     display_I = imadjust(ImageOnSphere);
%     greyColor = [.7 .7 .7];
%     [x,y,z]= sphere(30);
%     figure,
%     hold on
%     warp(x,y,z,(imcomplement(display_I)));
% %     warp(x,y,z,display_I);
%     surface(x,y,z,'FaceColor', 'none','EdgeColor',greyColor);
% %     surface(x,y,z,'FaceColor', display_I,'EdgeColor',greyColor);
% %     figure,
% %     hold on
% %     imagesc(ImageOnSphere);
% %     colorbar
% %     hold off
%     
    map = parula(256);
    minIval = min(min(ImageOnSphere));
    maxIval = max(max(ImageOnSphere));
    if(isPlotting)
        disp(['radius=', num2str(radius)]);
        disp(['min val=', num2str(minIval)]);
        disp(['max val=', num2str(maxIval)]);
    end
    
    display_I = 255 * (ImageOnSphere - minIval)/(maxIval - minIval);
%     greyColor = [.7 .7 .7];
%     colormap(map)
%     colorImg = ind2rgb(imcomplement(imadjust(ImageOnSphere)), map);
    colorImg = ind2rgb(uint8(display_I), map);
    
    ImageOnSphere = display_I /255; % is this better for accuracy
    
%     map = parula(256);
%     ImageOnSphere_rescaled = 255 * rescale(ImageOnSphere);
%     colorImg = ind2rgb(ImageOnSphere_rescaled, map);
    [x,y,z]= sphere(40);
    figure,
    hold on
    warp(x,y,z, colorImg, map );
    % surface(x,y,z,'FaceColor', 'none','EdgeColor',[0.7 0.7 0.7]);
    hold off
    axis equal;
    if(isRemoveAxesLabels)
        axis off
    else
        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
end
% %% Plot the vertices
% figure
% showPointCloud(desiredPointCloud);
% title('MagVals');
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% title('desiredPointCloud vertices at specific radius')
%
%
% figure
% showPointCloud(sphericalVertices);
% title('MagVals');
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% title('Spherical vertices at specific radius')

end