

%
% Rotate around x, y, z (counterclockwise when looking towards the origin)
%

function R = rotate_mat(x, y, z)
% It is assumed that this matrix is muplied from the hind, rotating objects
% along x-axis, y-axis, and z-axis in this order.

Rx = [      1       0       0; ...
            0  cos(x) -sin(x); ...
            0  sin(x)  cos(x)];

Ry = [ cos(y)       0  sin(y); ...
            0       1       0; ...
      -sin(y)       0  cos(y)];

Rz = [ cos(z) -sin(z)       0; ...
       sin(z)  cos(z)       0; ...
            0      0        1];

R = Rz*Ry*Rx;

return;