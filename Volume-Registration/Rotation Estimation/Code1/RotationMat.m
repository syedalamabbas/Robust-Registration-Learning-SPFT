
function [vertices, sph_verts, faces, fvec]=RotationMat(filename, confs)

load(filename);
[path,name,ext,ver] = fileparts(filename);

if ~exist('faces', 'var') || ~exist('vertices', 'var') || ~exist('sph_verts', 'var') || ~exist('fvec', 'var')
    disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
    return;
end
disp(sprintf('file name %s:',name)); 

degree = 1;

% rotate the parameter space
%disp('<< Rotate the parameter space >>');
[fvec, sph_verts] = param_rotate(fvec, vertices, sph_verts, faces, degree);


% rotate in the object space
% svs = [0 0 1; 1 0 0; 0 0 -1]; % north pole, intersection of dateline and equator, south pole
svs = [1 0 0; 0 -1 0; 0 0 1]; % x y z
Z = calculate_SPHARM_basis(svs, 1);

%disp('<< Rotate the object space >>');
vs = real(Z(:,2:4)*fvec(2:4,:));
R = object_rotate_R(vs);
[thetaX, thetaY, thetaZ] = factor_rot_xyz(R');
 disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',abs(thetaX/pi*180),abs(thetaY/pi*180),abs(thetaZ/pi*180)));

% X0=7;Y0=14;Z0=92;   % initial position
% % 
% X0=-59.81;     Y0=20.05 ;  Z0=166;
% XX=thetaX/pi*180-X0;
% YY=thetaY/pi*180-Y0;
% if ((XX> 12) || (YY > 12)||(XX< -12)||(YY< -12)) % thershold value 12
% ZZ=thetaZ/pi*180-Z0;
% else
% ZZ;
% end
%disp(sprintf('\n\nRotation ///xyz: X=%0.2f        Y=%0.2f       Z=%0.2f',XX,YY,ZZ));
% 
%  fvec = fvec*R'; 
%  vertices = vertices*R'; 



return;

%
% Factoring a Rotation Matrix as Rz*Ry*Rx (counterclockwise when looking towards the origin)
%

function [thetaX, thetaY, thetaZ] = factor_rot_xyz(R)

thetaY = asin(R(1,3));
if (thetaY < pi/2)
    if (thetaY > -pi/2)
        thetaX = atan2(-R(2,3),R(3,3));  % Rotation around z
        thetaZ = atan2(-R(1,2),R(1,1));
    else
        disp('WARNING (Factor Rotation): thetaY = -pi/2, not a unique solution, set thetaZ = 0');
        thetaX = -atan2(R(2,1),R(2,2));
        thetaZ = 0;
    end
else
    disp('WARNING (Factor Rotation): thetaY = pi/2, not a unique solution, set thetaZ = 0');
    thetaX = atan2(R(2,1),R(2,2));
    thetaZ = 0;
end

% disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',thetaX/pi*180,thetaY/pi*180,thetaZ/pi*180));

return;



%
% Parameter space rotation using degree 1 ellipsoid
% 

function [fvec, sph_verts,ZZ] = param_rotate(fvec, vertices, sph_verts, faces, degree)

% calculate matrix A
coeffs = fvec(2:4,:);
A(:,1) = (coeffs(1,:)-coeffs(3,:))';
A(:,2) = -(coeffs(1,:)+coeffs(3,:))'*i; 
A(:,3) = sqrt(2)*coeffs(2,:)';
A = real(A*sqrt(3)/(2*sqrt(2*pi)));


% SVD to find rotation and scaling matrics
%   [U,S,V] = svd(X) produces a diagonal matrix S and unitary matrices U and V so
%   that X = U*S*V'.
% need to rotate object space first.
[U,S,V] = svd(A);
%S the main axies that have half length
% extremum and saddle points
 % expts = U*S;

% set up rotation matrix
R = V(:,[3 2 1])';
   
if (det(R)<0)
    disp('WARNING (Parameter Space): Change back to pure rotation');
    R(2,:) = R(2,:)*(-1);
end

[thetaX, thetaY, thetaZ] = factor_rot_xyz(R');
%  [thetau]=DirectMethod(R');
% [thetau1]=EigenMethod(R');

% To multiply by -1
R = rotate_mat(thetaX, -thetaY, thetaZ);
%disp(sprintf('Rotation ///xyz: %0.2f %0.2f %0.2f',thetaX/pi,thetaY/pi,thetaZ/pi));

% the parameter space rotation
sph_verts = (R*sph_verts')';
% create new spharm descriptor (degree 1 is enough)
% [fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, 1, '', '');
% 
% % calculate [Estimating] the intersection point, north pole
% svs = [0 0 1; 1 0 0]; % north pole, intersection
% Z = calculate_SPHARM_basis(svs, 1);
% 
% vs = real(Z(:,2:4)*fvec(2:4,:));
% %disp(sprintf('vs: (%f, %f, %f)',vs));

% % north pole, intersection should be on the positive side 
% if vs(2,2)<0 || vs(1,3)<0
%     if vs(2,2)<0
%         Rfix = rotate_mat(0, 0, pi);
%     else
%         Rfix = eye(3);
%     end
%     if vs(1,3)<0
%         Rfix =  rotate_mat(pi, 0, 0)*Rfix;
%     end
% 	% the parameter space rotation
% 	sph_verts = (Rfix*sph_verts')';
%    	% create new spharm descriptor (degree 1 is enough)
%     
%    [fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, 1, '', '');
% end

 [fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, 1, '', '');

 return;


