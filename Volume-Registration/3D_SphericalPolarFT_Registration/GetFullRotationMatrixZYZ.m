function [ R ] = GetFullRotationMatrixZYZ( alpha, beta, gamma )
%GETFULLROTATIONMATRIX  computes the full 3 x 3 rotation matrix 
% given,
% alpha = rotation in radians about Z-axis  
% beta  = rotation in radians about Y-axis
% gamma = rotation in radians about Z-axis
Rz1 = [ cos(alpha)  -sin(alpha)   0
       sin(alpha)   cos(alpha)   0
           0             0       1 ];
Ry = [ cos(beta)     0        sin(beta)
         0           1            0
      -sin(beta)     0        cos(beta)];
Rz = [ cos(gamma)  -sin(gamma)   0
       sin(gamma)   cos(gamma)   0
           0             0       1 ];

R = Rz1 * Ry * Rz;            %% Full Rotation matrix based on Euler Angles
end

