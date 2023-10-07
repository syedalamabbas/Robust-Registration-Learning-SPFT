clear all;
close all;
clc;
load vertices.mat

data=vertices1;
origin=[0 0 0];
vxsize =[1 1 1];
 %function [vertices faces] = verticestosurface(data)
rangeX=floor(min(data(:,1))):1:ceil(max(data(:,1)));
rangeY=floor(min(data(:,2))):1:ceil(max(data(:,2)));
rangeZ=floor(min(data(:,3))):1:ceil(max(data(:,3)));

[X,Y]=meshgrid(rangeX,rangeY);

 X=round(imresize(X,0.21));
 Y=round(imresize(Y,0.21));
Z=round(griddata(data(:,1),data(:,2),data(:,3),X,Y,'cubic'));

% Z(any(isnan( Z),2),1:25) = 10;
mmin=min(min(Z));
 Z(:,:)=Z(:,:)-mmin;
 X(:,:)=X(:,:)-min(min(X));
 Y(:,:)=Y(:,:)-min(min(Y));
 
% ZZ=isnan(Z);
% Zmax=max(max(Z));
%   [M,N]=size(Z);
%  bim=zeros(M,N,255);
%  for k=1:Zmax
%                for ii=1:M
%                  for jj=1:N      
%                         if ((ZZ(ii,jj)== 1) || Z(ii,jj)<50)   %(Z(ii,jj)>=k))  %(ZZ(ii,jj)== 1)
% %                           X(ii+k,jj+k)=150;
% %                           Y(ii+k,jj+k)=120;
%                               Z(ii,jj)=50;
%                           bim(ii,jj,k)=1;
% %                        else
% %                            Z1(ii,jj)=0;
% %                            Z(ii,jj)=0;
%                         end
%                  end
%                end
% 
%  end         
  % surface(X,Y,Z)
      surf(X,Y,Z); 
%      
% %       [vertices,faces]=createsurface(bim);
%      
    [fvc3,fvc4] = u3d_pre; 
   vertices=fvc3.vertices;
   faces=fvc3.faces;
   
  % bim = simpvol( fvc4.vertices, fvc4.faces ); %one vector
 
%     bim=surf2vol(vertices,faces,rangeX,rangeY,rangeZ);
%    
%      figure, patch(isosurface(bim,0.1), 'Facecolor', [1 0 0]);
%      
 %________________________________________________________________________  
 %% Show the mesh
  figure, patch(fvc3,'FaceColor',[1 0 0]); axis square;

  % Convert the mesh to a voxelvolume
  Volume=polygon2voxel(fvc3,[50 50 50],'auto');
   Volume(:,:,1:10)=0;
  
  
  % Show x,y,z slices
  figure,
  subplot(1,3,1), imshow(squeeze(Volume(25,:,:)));
  title('Side view');
  subplot(1,3,2), imshow(squeeze(Volume(:,25,:)));
  title('Front view');
  subplot(1,3,3), imshow(squeeze(Volume(:,:,25)));
  title('Top view');
  %  Show iso surface of result
  figure, patch(isosurface(Volume,0.1), 'Facecolor', [1 0 0]);
 
   
% vertices(any(isnan( vertices),2),:) = [];
     %  vertices = vertices(faces(isfinite(faces)),:);
 PlotSurface(vertices,faces);
  %Calculate the smoothed version
  % FV2=smoothpatch(fvc,1,5);

   % Show the mesh and smoothed mesh
%    figure, 
%     subplot(1,2,1), patch(fvc,'FaceColor',[1 0 0],'EdgeAlpha',0);  view(3); camlight
%     subplot(1,2,2), patch(FV2,'FaceColor',[0 0 1],'EdgeAlpha',0); view(3); camlight
% 
%   
%  ZZ=isnan(vertices);
% rangeX=floor(min(vertices(:,1))):1:ceil(max(vertices(:,1)));
% rangeY=floor(min(vertices(:,2))):1:ceil(max(vertices(:,2)));
% rangeZ=floor(min(vertices(:,3))):1:ceil(max(vertices(:,3)));
% 
%  img=surf2vol(vertices,faces,rangeX,rangeY,rangeZ);
%  [vertices,faces]=createsurface(img);
% patch('Faces',faces,'Vertices',vertices);
%   %mesh_to_latex('C:\Users\Salah\Desktop\Eye gaze\matlab\tex\mesh',fvc.vertices,uint32(fvc.faces),fvc.facevertexcdata);