function PlotTwoRotationAngles( rotation_vector1, rotation_vector2 )
%PLOTTWOROTATIONANGLES Plots the rotation vectors as the spheres

figure,
[x,y,z]= sphere(30);
greyColor = [.7 .7 .7];
surface(x,y,z,'FaceColor', 'none','EdgeColor',greyColor);
axis equal;
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Plotting the two rotation vectors to determine if they are close or not')
view(-19,26)

LinePoints = -1.4:.02:1.4;                % Points around and covering the unit sphere
pointCloudLineX = [LinePoints', zeros(size(LinePoints))',zeros(size(LinePoints))'];
pointCloudLineY = [zeros(size(LinePoints))',LinePoints',zeros(size(LinePoints))'];
pointCloudLineZ = [zeros(size(LinePoints))',zeros(size(LinePoints))',LinePoints'];
pointCloudLine = [pointCloudLineX ;pointCloudLineY; pointCloudLineZ];
ptCloud = pointCloud(pointCloudLine);

t = eye(4);
t(1:3,1:3)= GetFullRotationMatrixZYZ( rotation_vector1(1), rotation_vector1(2), rotation_vector1(3));
tformRotated = affine3d(t);
ptCloudTransformed =  pctransform(ptCloud,tformRotated);
pointCloudLineTransformed = ptCloudTransformed.Location;
PlotMarkerLinesAtOrigin(pointCloudLineTransformed, 36, [0 1 0], 'o');


t = eye(4);
t(1:3,1:3)= GetFullRotationMatrixZYZ( rotation_vector2(1), rotation_vector2(2), rotation_vector2(3));
tformRotated = affine3d(t);
ptCloudTransformed =  pctransform(ptCloud,tformRotated);
pointCloudLineTransformed = ptCloudTransformed.Location;
PlotMarkerLinesAtOrigin(pointCloudLineTransformed, 36,[1 0 0],'s');     % Plotting marker lines over the images obtained from the above function


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

