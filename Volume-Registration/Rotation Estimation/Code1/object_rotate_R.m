%
% rotation matrix in object space
%

function R = object_rotate_R(vs)

% fix north pole
[PHI,THETA] = cart2sph(vs(1,1),vs(1,2),vs(1,3));
% if phi < 0 shift it by 360
ind = find(PHI<0); PHI(ind) = PHI(ind)+2*pi;
%  shift THETA by + 90
THETA = pi/2-THETA;
alpha = -PHI; beta = -THETA;
R = rotate_mat(0, beta, 0)*rotate_mat(0, 0, alpha);
% Rotate the the object space 
vs = vs*R';
% fix intersection;
[PHI,THETA] = cart2sph(vs(2,1),vs(2,2),vs(2,3));
gamma = -PHI;
R1 = rotate_mat(0, 0, gamma); R = R1*R;
vs = vs*R1';

return;