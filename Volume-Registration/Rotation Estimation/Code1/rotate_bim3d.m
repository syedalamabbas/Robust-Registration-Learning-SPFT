
function [bim, flag,theta_r] = rotate_bim3d(inNames,theta_r)

  load(inNames{1});
    if exist('bim','var')
        blob1 = bim;
    end

 DIM = size(blob1);
x=40;  % more than 55 out of memory error

% make a work area to rotate a round Z
d = DIM+x; blob = zeros(d);
blob((x/2+1):d(1)-x/2,(x/2+1):d(2)-x/2,(x/2+1):d(3)-x/2) = blob1;

blob_center = (size(blob) + 1) / 2;
%blob_center = [ 25 28 15];


T1 = [1 0 0 0; 0 1 0 0;  0 0 1 0; -blob_center 1];

% 
% theta_r.x= input('Enter the rotation angle around X: ');
   thetaxr=theta_r.x*(pi/180);
%   theta_r.x=(rand_int(-177,179,1));
%    thetaxr=theta_r.x*(pi/180);

if (thetaxr > 1.5708)
    thetax=thetaxr - 1.5708;
    flag.x=-1;
elseif (thetaxr < -1.5708)
    thetax=thetaxr + 1.5708;
    flag.x=1;
else
    thetax=thetaxr;
    flag.x=0;
end 

%   
%   theta_r.y= input('Enter the rotation angle around Y: ');
    thetayr=-theta_r.y*(pi/180);
%  theta_r.y=(rand_int(-177,179,1));
%  thetayr=-theta_r.y*(pi/180);
 if (thetayr > 1.5708)
    thetay=thetayr - 1.5708;
    flag.y=1;
elseif (thetayr < -1.5708)
    thetay=thetayr + 1.5708;
    flag.y=-1;
 else
    thetay=thetayr;
    flag.y=0;
 end 


%   theta_r.z= input('Enter the rotation angle around Z: ');
   thetazr=theta_r.z*(pi/180);
%  theta_r.z=(rand_int(-177,179,1));
%  thetazr=-theta_r.z*(pi/180);
%  
if (thetazr > 1.5708)
    thetaz=thetazr - 1.5708;
    flag.z=-1;
elseif (thetazr < -1.5708)
    thetaz=thetazr + 1.5708;
    flag.z=1;
else
   thetaz=thetazr;
    flag.z=0;
end 


T2y = [cos(thetay)  0      -sin(thetay)   0;   0  1 0   0;  sin(thetay)   0   cos(thetay)   0;  0   0   0   1]; %Round Y
T2z = [cos(thetaz)   -sin(thetaz) 0  0; sin(thetaz)   cos(thetaz)  0  0;0  0   1   0;  0   0   0   1]; %
T2x = [1  0  0  0; 0  cos(thetax)   -sin(thetax)  0; 0 sin(thetax)   cos(thetax)   0; 0   0   0   1];  % Around X
T2=T2z * T2y * T2x;    %ZYX
T3 = [1  0  0  0; 0 1  0  0; 0  0  1  0; blob_center 1];
%T3 = [1  0  0  0; 0 1  0  0; 0  0  1  0; 0  0  0  1];
T = T1 * T2 * T3;

tform = maketform('affine', T);
% The tform struct should map the blob center to itself.
tformfwd(blob_center, tform);
% --------------------------------------------------------------------------
% The type of interpolation 
R = makeresampler('linear', 'fill');
% TDIMS_A & TDIMS_B specifies how the dimensions of the output array 
%  correspond to the dimensions of the spatial transformation. 
TDIMS_A = [1 2 3];
TDIMS_B = [1 2 3];
% TSIZE_B is the size of the output array.
 TSIZE_B = size(blob);
  %TSIZE_B = [120 120 255];
% TMAP_B is unused when you have a tform struct. Just specify it to be empty.
   TMAP_B = [];
% F specifies the values to use outside the boundaries of the input array.
   F = 0;
% Call tformarray to transform the blob
bim = tformarray(blob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);