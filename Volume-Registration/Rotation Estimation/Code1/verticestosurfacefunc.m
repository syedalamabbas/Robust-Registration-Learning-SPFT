

function [vertices faces] = verticestosurfacefunc(data)
rangeX=floor(min(data(:,1))):1:ceil(max(data(:,1)));
rangeY=floor(min(data(:,2))):1:ceil(max(data(:,2)));
%rangeZ=floor(min(data(:,2))):1:ceil(max(data(:,2)));

[X,Y]=meshgrid(rangeX,rangeY);
Z=round(griddata(data(:,1),data(:,2),data(:,3),X,Y,'v4'));
mmin=min(min(Z));
Z(:,:)=Z(:,:)-mmin;

 X=round(imresize(X,0.5));
 Y=round(imresize(Y,0.5));
 Z=round(imresize(Z,0.5));
 mmin=min(min(Z));
ZZ=isnan(Z);
% Zmax=max(max(Z));
[M,N]=size(Z);
%bm=zeros(M,N,256);
              for ii=1:M
                for jj=1:N      
                       if (ZZ(ii,jj)== 1)
                           Z(ii,jj)=mmin;
                           X(ii,jj)=mmin;
                           Y(ii,jj)=mmin;
                       
                       end
                end
              end

          
   %surface(X,Y,Z)
     surf(X,Y,Z); 
  fvc = u3d_pre; 
 vertices=fvc.vertices;
 faces=fvc.faces;
%   %mesh_to_latex('C:\Users\Salah\Desktop\Eye gaze\matlab\tex\mesh',fvc.vertices,uint32(fvc.faces),fvc.facevertexcdata);