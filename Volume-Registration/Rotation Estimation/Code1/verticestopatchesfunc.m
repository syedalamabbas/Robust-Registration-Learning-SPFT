

function [vertices faces] = verticestopatchesfunc(data)
rangeX=floor(min(data(:,1))):1:ceil(max(data(:,1)));
rangeY=floor(min(data(:,2))):1:ceil(max(data(:,2)));
%rangeZ=floor(min(data(:,2))):1:ceil(max(data(:,2)));

[X,Y]=meshgrid(rangeX,rangeY);
 
 X=round(imresize(X,0.5));
 Y=round(imresize(Y,0.5));
%  Z=round(imresize(Z,0.5));
Z=round(griddata(data(:,1),data(:,2),data(:,3),X,Y,'cubic'));

mmin=min(min(Z));
 Z(:,:)=Z(:,:)-mmin;
 X(:,:)=X(:,:)-min(min(X));
 Y(:,:)=Y(:,:)-min(min(Y));
%

Zmax=max(max(Z));
 [M,N]=size(Z);

  bim=zeros(M,N,255);
 for k=Zmax-50:Zmax-40
              for ii=1:M
                for jj=1:N      
                       if ((Z(ii,jj)>=k))  %(ZZ(ii,jj)== 1)||
%                           X(ii+k,jj+k)=150;
%                           Y(ii+k,jj+k)=120;
%                           Z(ii,jj)=80;
                          bim(ii,jj,k)=1;
%                        else
%                            Z1(ii,jj)=0;
%                            Z(ii,jj)=0;
                       end
                end
              end

 end  
 vdata=data;
 fdata=0;
   
     [vertices,faces]=createsurface(bim,vdata,fdata);          
   %surface(X,Y,Z)
%      surf(X,Y,Z); 
%   fvc = u3d_pre; 
%  vertices=fvc.vertices;
%  faces=fvc.faces;
%   %mesh_to_latex('C:\Users\Salah\Desktop\Eye gaze\matlab\tex\mesh',fvc.vertices,uint32(fvc.faces),fvc.facevertexcdata);