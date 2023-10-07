
function [roll pitch yaw]=calc_pose(inNames)

%% (4) Estimate the Rotation Matrix

numSbj = length(inNames);


h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file = inNames{i};

            load(file);
            [path,name,ext,ver] = fileparts(file);

             if ~exist('faces', 'var') || ~exist('vertices', 'var') || ~exist('sph_verts', 'var') || ~exist('fvec', 'var')
             disp('Something are missing from faces, vertices, spherical vertices');
             return;
             end
    disp(sprintf('file name %s:',name)); 
    
%% rotate the parameter space
%disp('<< Rotate the parameter space >>');

% calculate the matrix A from frequency coefficients
coeffs = fvec(2:4,:);
A(:,1) = (coeffs(1,:)-coeffs(3,:))';
A(:,2) = -(coeffs(1,:)+coeffs(3,:))'*sqrt(-1); 
A(:,3) = sqrt(2)*coeffs(2,:)';
A = real(A*sqrt(3)/(2*sqrt(2*pi)));

% SVD to decompose the Matrix A
[U,S,V] = svd(A);
% set up rotation matrix
R = V(:,[3 2 1])';
%Calculate the rotation angles
[thetaX, thetaY, thetaZ] = factor_rot_xyz(R');

% To multiply by -1
R = rotate_mat(thetaX, -thetaY, thetaZ);

% the parameter space rotation
sph_verts = (R*sph_verts')';
%Create the ellipose
 [fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, 1, '', '');

%% rotate  the object space
%disp('<< Rotate the object space >>');
p = [1 0 0; 0 -1 0; 0 0 1]; % x y z
Z = calculate_SPHARM_basis(p, 1);

vs = real(Z(:,2:4)*fvec(2:4,:));


%  Rotate the object space and set up rotation matrix
R = object_rotate_R(vs);
%Calculate the rotation angles
[thetaX, thetaY, thetaZ] = factor_rot_xyz(R');

% Shift the Z 0.7
% if ((abs(thetaZ) < 83*pi/180))
%      thetaZ=abs(thetaZ)+8*pi/180;  
% else
%          thetaZ=abs(thetaZ)-8*pi/180; 
% end 
% 
 if (abs(thetaZ) > pi/2)
     thetaZ=pi-abs(thetaZ); 
 end
  thetaZ=abs(thetaZ)-( 40*pi/180);  
  if (abs(thetaX) > pi/2)
     thetaX=pi-abs(thetaX); 
 end
%  if (thetaX < 0)
%     thetaX=thetaX+2*pi;  
% else
%     thetaX=thetaX-2*pi;
% end   
% disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',abs(thetaX/pi*180),abs(thetaY/pi*180),abs(thetaZ/pi*180)));
  disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',round(thetaX/pi*180),round(thetaY/pi*180),round(thetaZ/pi*180)+0));

 clear('vertices', 'sph_verts', 'faces', 'fvec');
    waitbar(i/numSbj)
 end
 close(h);
roll=thetaX; pitch=thetaY;  yaw=thetaZ;