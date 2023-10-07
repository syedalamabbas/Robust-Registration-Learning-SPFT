% Create surfaces from 3-D voxel data  

% Input 
% voxels(nrow,ncol,nslice)  Ids of all voxels : nrow=ny ncol=nx nslice=nz
% sf(1,3)                   Voxel scale factors for x,y,z dimensions
%                           sf=[1 1 1] if problem does not have any dimensional scale
% offset(1,3)               Offset of crd of center of voxel(1,1,1) from global origin
%                           offset=[0 0 0] if global origin is the center of voxel(1,1,1)
% id                        Voxel id whose boundary is to be determined

% Output
% faces(nfaces,3)           Connectivity of boundary faces : vertices are numbered contigously from 1 to nvtx
% vertex(nvtx,dim)          Crd of boundary vertices : dim=3

% See fnc=voxel_bnd_faces for details

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [fighndl]= PlotSurfaceT(vertices,faces,filename,theta,theta_r)
id=1;
% Colors for plotting
red=[1.0 0.0 0.0];  
green=[0.5 1.0 0.5];  
cyan=[0.0 1.0 1.0];  
blue=[0.0 0.0 1.0]; 
yellow=[1.0 1.0 0.0];  
magenta=[1.0 0.0 1.0];
orange=[1.0 0.5 0.0];  
brown=[1.0 0.5 0.5]; 
purple=[0.5 0.5 1.0]; 
black=[0.0 0.0 0.0]; 
noclr='none';


% 
%  load matf_imagesnifti  % image file
%  voxels=Images{1};  % voxels
% load imagebim133_55t  % image file
% voxels=bim; %imresize(bim,0.5);  %bim; %{1};  % voxels
% 
% voxid=1;  % voxel domain ids

% Scale factors for scaling each dimension
% sf=[0.5508 0.5508 0.6000];
sf=[1 1 1];
% 1 voxel along x axis=0.5508 units of length
% 1 voxel along y axis=0.5508 units of length
% 1 voxel along z axis=0.6 units of length
  
% Offset of crd of center of voxel(1,1,1) from global origin
offset=[0 0 0];  % center of voxel(1,1,1) = global origin

facecolor=[  % color of each domain
red
green
cyan
blue
yellow
]; 

% edgecolor='none';
edgecolor=[0 0 0];  % black

na=[ 'Object: ' filename  ]; %'   Thetha_x=' num2str(theta.x) '   Theta_y='  num2str(theta.y) '  Theta_z='  num2str(theta.z)];
tit=({['Object: ' filename '   Thetha_x=' num2str(theta.x) '   Theta_y='  num2str(theta.y) '  Theta_z='  num2str(theta.z)],[ ] , ['      Ground Truth       Thetha_x=' num2str(theta_r.x) '   Theta_y='  num2str(theta_r.y) '  Theta_z='  num2str(theta_r.z)]});
fighndl=[];

% -------------------------------------------------------------------------
% Create surfaces for each voxel domain
% In this example faces and vertex are overwritten
% Store them in a cell array or structure if it is required to process them

% for id=voxid
%     id
    

    % Plot surface
    [fighndl]=plotsurf(faces,vertices,facecolor(id,:),edgecolor,na,tit,fighndl);  
%     name= ['..\3Ddataset\in_images\' filename '.png'];
%       saveas (fighndl ,name ); 
     % fighndl=[];
% end
