
function []=make_mesh_structure(MeshZ)
% clear all;
% close all;
% clc

%I = imread('depthmap.png'); %read in a test image
%I = imnoise(I, 'salt & pepper');
%figure(1); imshow(I, [])
%[rows,cols] = size(I);
%[MeshX, MeshY]=meshgrid(1:cols,1:rows);
%MeshZ=rgb2gray(I);
MeshZ=double(MeshZ); 
[rows,cols]=size(MeshZ);
%MeshZ=(max(max(MeshZ))-min(min(MeshZ)))/min(min(MeshZ));
MeshZ= (MeshZ/max(max(MeshZ)))*30;
for i=1:rows
      for j=1:cols
            if MeshZ(i,j)<=10;
                  MeshZ(i,j)=NaN;
            end
      end
 end
 
if (sum(sum(~isnan(MeshZ))) > 20)
    X = 1:rows;
    Y = (1:cols)';
    Z = MeshZ;
        fvc3 = surf2patch(X,Y,Z, 'triangles' );
       %  fvc4 = surf2patch(X,Y,Z,realcolor);
else
    fvc3.vertices            = [];
    fvc3.faces               = [];
    fvc3.facevertexcdata     = [];
    disp('No surface graph found.');
     fvc4.vertices            = [];
    fvc4.faces               = [];
    fvc4.facevertexcdata     = [];
    disp('No surface graph found.');

end;
%mesh(MeshZ);
%  fvc = u3d_pre; 
 vertices=fvc3.vertices;
 faces=fvc3.faces;
  figure, drawMesh(vertices,faces,'wire');
   filename= 'filename';
     theta.x=0;
     theta.y=0;
     theta.z=0;
     %filename= 'name ';
     PlotSurface(vertices,faces,filename,theta);
     % patch('Faces',faces,'Vertices',vertices);
      % Convert the mesh to a voxelvolume
 
     % Volume=polygon2voxel(fvc3,[50 50 50],'auto');
       [bim] = verticestovolumefunc(vertices,faces,filename); % save new bim
       
%% (1) Input data directory
Spharm3DModelDir = '../';
CodeDir = [Spharm3DModelDir 'code1'];
addpath(CodeDir);

DataFolder = '3Ddataset';
InDataDir = [Spharm3DModelDir DataFolder '/dataTest'];


%% (2) List input file names


    inFiles = dir([InDataDir '/*.mat']); inNames={};
     for ii=1:length(inFiles)
          inNames{end+1} = [InDataDir '/' inFiles(ii).name];
     end 
%% (3) Display input objects (optional)
% 
    %Available values for Space- 'object';'param';'both'
dispConfs.Space = 'object';
    % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    %    'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'orig';  
    % Available values for Shape- 'solid';'mesh';'both'
dispConfs.Shade = 'both';
    % Available values for Overlay- 'none';'adc_paramap'
dispConfs.Overlay = 'none';
    % Available values for Export- 'screen';'png';'both'
dispConfs.Export = 'both';
dispConfs.Degree = [];
dispConfs.Template = '';
% the name of the orginal image to be stored
dispConfs.name=filename;

SpharmMatUtilDisplayObjsPS(dispConfs, inNames, CodeDir);

%%
  
%V=[MeshX(:), MeshY(:), MeshZ(:)];

% create and interpolate the mesh to be 10x10
%[V2, F2]=create_structured_mesh_sorted(V, 100);

% display both meshes
%figure, drawMesh(V1,F1,'wire');
%figure, drawMesh(V2,F2,'wire');

end

