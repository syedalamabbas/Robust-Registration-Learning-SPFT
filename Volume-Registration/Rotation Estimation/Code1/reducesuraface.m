clear all;
close all;
clc;
load vertices.mat

data=vertices1;

 %function [vertices faces] = verticestosurface(data)
rangeX=floor(min(data(:,1))):1:ceil(max(data(:,1)));
rangeY=floor(min(data(:,2))):1:ceil(max(data(:,2)));
%rangeZ=floor(min(data(:,3))):1:ceil(max(data(:,3)));

[X,Y]=meshgrid(rangeX,rangeY);
% 
 X=round(imresize(X,0.05));
 Y=round(imresize(Y,0.05));
Z=round(griddata(data(:,1),data(:,2),data(:,3),X,Y)); %,'cubic'));
mmin=min(min(Z));
 Z(:,:)=Z(:,:)-mmin;
 X(:,:)=X(:,:)-min(min(X));
 Y(:,:)=Y(:,:)-min(min(Y));
 
% ZZ=isnan(Z);
Zmax=max(max(Z));
 [M,N]=size(Z);
 bm=zeros(M,N,(Zmax-50):Zmax);
 for k=(Zmax-50):Zmax
              for ii=1:M
                for jj=1:N      
                       if ((Z(ii,jj)>=k))  %(ZZ(ii,jj)== 1)||
%                           X(ii+k,jj+k)=150;
%                           Y(ii+k,jj+k)=120;
%                           Z(ii,jj)=80;
                          bm(ii,jj,k)=1;
%                        else
%                            Z1(ii,jj)=0;
%                            Z(ii,jj)=0;
                       end
                end
              end

 end         
  % surface(X,Y,Z)
     surf(X,Y,Z); 
     
     [vertices,faces]=createsurface(bm);
     
%   fvc = u3d_pre; 
%  vertices=fvc.vertices;
%  faces=fvc.faces;
%  ZZ=isnan(vertices);
% rangeX=floor(min(vertices(:,1))):1:ceil(max(vertices(:,1)));
% rangeY=floor(min(vertices(:,2))):1:ceil(max(vertices(:,2)));
% rangeZ=floor(min(vertices(:,3))):1:ceil(max(vertices(:,3)));
% 
% img=surf2vol(vertices,faces,rangeX,rangeY,rangeZ);
%patch('Faces',faces1,'Vertices',vertices1);
%   %mesh_to_latex('C:\Users\Salah\Desktop\Eye gaze\matlab\tex\mesh',fvc.vertices,uint32(fvc.faces),fvc.facevertexcdata);