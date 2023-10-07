close all


%% My solution
B = 8; % Some Even number of bandwidth
l=0:2*B-1;                         % 2-Sphere of 2B x 2B equiangular grid 
k = 0:2*B-1;
theta = (2*l+1)*180/(4*B);
phi = 2*180*k/(2*B);

r = 1;
greyColor = [.7 .7 .7];

figure,
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
axis tight % ([-1.1 1.1 -1.1 1.1 -1.1 1.1])
view(-60,20);
grid on
hold on
for l = 1:2*B
    for k = 1:2*B
        X = r .* sind(theta(l)) .* cosd(phi(k));                     %% Note the conventions are different that is used  here to get the equiangular grid
        Y = r .* sind(theta(l)) .* sind(phi(k));
        Z = r .* cosd(theta(l)); 
        s = scatter3(X,Y,Z,'o','LineWidth',1.2);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = 'b';
        drawnow
        if(l == 1)
            plot3( r .* sind([theta 180 ]) .*cosd(phi(k)), r .* sind([theta 180 ]) .* sind(phi(k)), r .* cosd([theta 180 ]),'color', greyColor);
        end
    end
    plot3( r .* sind(theta(l)) .*cosd([phi 360 ]), r .* sind(theta(l)) .* sind([phi 360 ]), r .* cosd(theta(l)).* ones(size([phi 360 ])),'color', greyColor);  % Fix theta , all phi's
end

r = 2;
for l = 1:2*B
    for k = 1:2*B
        X = r .* sind(theta(l)) .* cosd(phi(k));                     %% Note the conventions are different that is used  here to get the equiangular grid
        Y = r .* sind(theta(l)) .* sind(phi(k));
        Z = r .* cosd(theta(l)); 
        s = scatter3(X,Y,Z,'o','LineWidth',1.2);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = 'r';
        drawnow
        if(l == 1)
            plot3( r .* sind([theta 180 ]) .*cosd(phi(k)), r .* sind([theta 180 ]) .* sind(phi(k)), r .* cosd([theta 180 ]),'color', greyColor);
        end
    end
    plot3( r .* sind(theta(l)) .*cosd([phi 360 ]), r .* sind(theta(l)) .* sind([phi 360 ]), r .* cosd(theta(l)).* ones(size([phi 360 ])),'color', greyColor);  % Fix theta , all phi's
end

r = 3;
for l = 1:2*B
    for k = 1:2*B
        X = r .* sind(theta(l)) .* cosd(phi(k));                     %% Note the conventions are different that is used  here to get the equiangular grid
        Y = r .* sind(theta(l)) .* sind(phi(k));
        Z = r .* cosd(theta(l)); 
        s = scatter3(X,Y,Z,'o','LineWidth',1.2);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = 'g';
        drawnow
        if(l == 1)
            plot3( r .* sind([theta 180 ]) .*cosd(phi(k)), r .* sind([theta 180 ]) .* sind(phi(k)), r .* cosd([theta 180 ]),'color', greyColor);
        end
    end
    plot3( r .* sind(theta(l)) .*cosd([phi 360 ]), r .* sind(theta(l)) .* sind([phi 360 ]), r .* cosd(theta(l)).* ones(size([phi 360 ])),'color', greyColor);  % Fix theta , all phi's
end
hold off

% print SphericalPolarGrid_123 -dpdf

%% Plotting multilayered spherical images
figure,
% axis equal

xlabel('x')
ylabel('y')
zlabel('z')
axis tight % ([-1.1 1.1 -1.1 1.1 -1.1 1.1])
view(-60,30);
grid on
hold on
r = 1;
[meshX,meshY,meshZ] =meshgrid(1:2*B,1:2*B,r);
surf(meshX,meshY,meshZ);
alpha .3
for l = 1:2*B
    for k = 1:2*B
        X = l;
        Y = k;
        Z = r; 
        s = scatter3(X,Y,Z,'o','LineWidth',1.2);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = 'b';
        drawnow
    end
end
r = 2;
[meshX,meshY,meshZ] =meshgrid(1:2*B,1:2*B,r);
surf(meshX,meshY,meshZ);
alpha .3
for l = 1:2*B
    for k = 1:2*B
        X = l;
        Y = k;
        Z = r; 
        s = scatter3(X,Y,Z,'o','LineWidth',1.2);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = 'r';
        drawnow
    end
end
r = 3;
[meshX,meshY,meshZ] =meshgrid(1:2*B,1:2*B,r);
surf(meshX,meshY,meshZ);
alpha .3
for l = 1:2*B
    for k = 1:2*B
        X = l;
        Y = k;
        Z = r; 
        s = scatter3(X,Y,Z,'o','LineWidth',1.2);
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = 'g';
        drawnow
    end
end
hold off
% print MultilayerSphericalImages_123 -dpdf
