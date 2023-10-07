clear all;
close all;
clc;
load vertices.mat

data=vertices1;

[vertices faces] = verticestopatchesfunc(data);

maxfn = 10^4; % Maximum number of faces
% adjust face numbers
fn = size(faces,1);
% if fn>maxfn
    disp(sprintf('Reduce face number from %d to %d',fn,maxfn));
    [faces,vertices] = reducepatch(faces,vertices,maxfn);
% end
PlotSurface(vertices,faces);
%[vertices,faces]=createsurface(bim,vertices,faces);  
%[vertices faces] = verticestopatchesfunc(vertices);
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