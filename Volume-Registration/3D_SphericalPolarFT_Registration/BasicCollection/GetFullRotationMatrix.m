function [ R ] = GetFullRotationMatrix( alpha, beta, gamma )
%GETFULLROTATIONMATRIX  computes the full 3 x 3 rotation matrix 
% given,
% alpha = rotation in radians about X-axis  
% beta  = rotation in radians about Y-axis
% gamma = rotation in radians about Z-axis
Rx = [ 1           0              0
       0       cos(alpha)  -sin(alpha)
       0       sin(alpha)   cos(alpha)];
Ry = [ cos(beta)    0      -sin(beta)
         0          1              0
      sin(beta)     0       cos(beta)];
Rz = [ cos(gamma)  -sin(gamma)   0
       sin(gamma)   cos(gamma)   0
           0             0       1 ];

R = Rx * Ry * Rz;            %% Full Rotation matrix based on Euler Angles
end

