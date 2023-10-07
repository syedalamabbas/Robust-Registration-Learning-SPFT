%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D:\Eye gaze\Spharm3DModel\3Ddataset\dataTest\test  XYZ positive -xy-z
%negative objectxyz  object-xy-z
clc
%% you have to create the average object first as Template
%%
clear all;
close hidden all;
clear hidden all;
count=0;
InDataDir = '..\';


 DataFolder = 'benchmark\3DObject_Dataset\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];
  deg1=15;
% Flower m1000 m1003  m1001 m1019 m1027 
% human  m119 m121 m122 m124 m156 194 184 163
% head 343 352 350 356 341
% human1  233 241 243 245 246 240 236
% bottel 482 485 486 490 492 
% fish 55 56 57 60  81
% gaithar 628 634 636 637 639 625
% chair 812 813 814 816 820  819   10
% air-plane m1194  1162  1161  1157  1165
% plum tree  m1089  m1091  m1096 m1094
% tree   m1046  1053  1054  1048 1040
% dbbd 1412 1413  1415 1418 1421 1422
% plane 1252 1253 1261  1262 1263
% ship 1430 1432  1433 1443 1427 1436  1446
% bike 1473 1475 1476
% mtorcyle  1481 1483 1484
% race car 1502 1510 1512 1514  1515
% car 1517 1518 1530 1531 1521 1525  26 27 33

jj=[1000 1003  1001 1019 1027 ];
%jj(:)=jj(:)+1;
  %for count=1:5  %40 19-24car
 %inFiles = dir([InDataDir   DataFolder   'm'  jj '.ply']);

   
     inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
     for ii=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
     end

    %*******************************************************************  
%% (1) Input data directory     
 %for i=1:length(inFiles)     
for i=1:5
     load(inNames_obj{i});

      [pa,name,ex]=fileparts(inNames_obj{i});     
      filetemplate=name;
     
      vertices= verticesOrg;
      faces=facesOrg;
      
      
   origin=[0 0 0];
   vxsize =[1 1 1];
confs.Connectivity='(6+,18)'; 
confs.Epsilon = 1.5000;
confs.OutDirectory = [InDataDir DataFolder '/dataFix'];
confs.name=name;
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

[bim] = verticestovolumefunc(vertices,faces);

[bim, new_name] = fix_holes(bim, confs);

[vertices, faces] =  gen_surf_data(bim,origin,vxsize);
    % max_degree0 = sqrt(size(vertices,1))-1
    %% % (4) Make 3D Rotation

     [vertices,blob_center] = center(vertices);
   % [vertices] = rotate_obj3d(vertices,theta_r);

%   temp=vertices(:,1);
%   vertices(:,1)=vertices(:,2);
%   vertices(:,2)=temp;
%%  
 PlotSurface1(vertices,faces);
%save(file,'vertices','faces');
    % close all;


%% (4.1) Perform Surface Meshes (CALD)

%Create spherical parameterization for surface meshes  **Mapping**
%// Read Oject
confs.MeshGridSize = 50;
confs.MaxSPHARMDegree =deg1;  % 6
confs.Tolerance = 2;
confs.Smoothing = 2;
confs.Iteration = 10; %100;
confs.LocalIteration = 1; %10;
    % Available values for t_major- 'x';'y'
confs.t_major = 'x';
    % Available values for SelectDiagonal- 'ShortDiag';'LongDiag'
confs.SelectDiagonal = 'ShortDiag';

confs.OutDirectory = [InDataDir DataFolder '\dataSurface'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

[sph_verts, new_name] = initParamCALD(vertices, faces, new_name, confs);

%///////////////////////////////inNames_obj   no rotation
%% (5) Perform SPHARM-MAT Expansion
% Calculate SPHARM coefficients
confs.MaxSPHARMDegree = deg1;   % 15
confs.OutDirectory = [InDataDir DataFolder '/dataExp'];
%///////////////////////////////////////////////////////////////////////
%  vertnum = size(sph_verts,1);
%  maxDeg= max(1, floor(sqrt(vertnum)*1/2));
% confs.MaxSPHARMDegree = maxDeg;
maxDeg=deg1;
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
outNames = {};

 %outNames = SpharmMatExpansion(confs, outNames, 'ExpLSF');
 [fvec, deg, Z, outNames{end+1}]=create_SPHARM_des_LSF([],[],[],maxDeg,new_name,deblank(char(confs.OutDirectory)));


%% (4) Perform FOE Alignment
    % Available values for CPoint- 'x';'y';'z'
% confs.CPoint = 'y';
%     % Available values for NPole- 'x';'y';'z'
% confs.NPole = 'z';
% confs.MaxSPHARMDegree = deg1;
% confs.OutDirectory = [InDataDir DataFolder 'data_reg'];
% 
% if ~exist(confs.OutDirectory,'dir')
%     mkdir(confs.OutDirectory);
% end
% 
% outNames = SpharmMatAlignment(confs, outNames, 'AligFOE');   

%% (7)estimate Rotation

confs.Template = [InDataDir DataFolder '/dataOutput'  '/atlas.mat'];
confs.MaxSPHARMDegree = deg1;   %12   depends on the vertices number deg = max(1, floor(sqrt(vertnum)*1/2));
confs.GroupAlpha = 100;
    % Available values for NormalizeSize- 'No';'Yes'
confs.NormalizeSize = 'No';
confs.BaseRes = 1;
confs.HierarchyStep = 1;
confs.HierarchyDepth = 3;
confs.Top_K = 1;
confs.GammaRes = 2;
confs.thetha(1,1:3)=0;
confs.OutDirectory = [InDataDir  DataFolder 'OutDataPose'];
confs.theta.x=0;
confs.theta.y=0;
confs.theta.z=0;

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
method='Eigenvector';
%ICP_Eigenvector, Eigenvector, FOE
[rmsd, theta] = SpharmMethod(confs, outNames, method);  
data(k,1)=rmsd;
data(k,2)=theta.x;
data(k,3)=theta.y;
data(k,4)=theta.z;
 
     end
%        end %IF STATMENT
%    end
 disp(' Average Objects')

%% (1) Input data directory
inFiles = dir([InDataDir  DataFolder 'data_reg/alignParam'  '/*_prm.mat']); inNames={};


%% (2) List input file names
for i=1:length(inFiles)
    inNames{end+1} = [InDataDir DataFolder 'data_reg/alignParam/'  inFiles(i).name];
end

 
%% (4) Average objects
confs.OutputName='atlas.mat'; 
confs.OutDirectory = [InDataDir    DataFolder 'dataOutput'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

SpharmMatUtilAverageObjs(confs, inNames);

 %end

%% (5) Display output objects (optional)
    %Available values for Space- 'object';'param';'both'
dispConfs.Space = 'object';
    % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    %    'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'icosa4';  
    % Available values for Shape- 'solid';'mesh';'both'
dispConfs.Shade = 'both';
    % Available values for Overlay- 'none';'adc_paramap'
dispConfs.Overlay = 'adc_paramap';
    % Available values for Export- 'screen';'png';'both'
dispConfs.Export = 'png';
dispConfs.Degree = [];
dispConfs.Template = '';

inFiles = dir([confs.OutDirectory '/*_als.mat']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);

TinNames={};
TinNames{1} = [confs.OutDirectory '/atlas.mat']; 

SpharmMatUtilDisplayObjs(dispConfs, TinNames, CodeDir);

%% (7) Estimate the rotation

% DataFolder = '3Ddataset\datatest\test';
 confs.Template = [confs.OutDirectory '/atlas.mat'];

 
confs.MaxSPHARMDegree = deg1;   %12   depends on the vertices number deg = max(1, floor(sqrt(vertnum)*1/2));
confs.GroupAlpha = 100;
    % Available values for NormalizeSize- 'No';'Yes'
confs.NormalizeSize = 'No';
confs.BaseRes = 1;
confs.HierarchyStep = 1;
confs.HierarchyDepth = 3;
confs.Top_K = 1;
confs.GammaRes = 2;
confs.thetha(1,1:3)=0;
%confs.OutDirectory = [InDataDir  DataFolder '/OutDataPose'];
confs.theta.x=0;
confs.theta.y=0;
confs.theta.z=0;

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
method='Eigenvector';

%ICP_Eigenvector, Eigenvector, FOE
[ rmsd, theta] = SpharmMethod(confs, inNames, method);   
%[rmsd] = SpharmMatAlignment(confs, outNames, 'AligSHREC');
%[ rmsd, confs] = SpharmMat(confs, outNames, 'Eigenvector');   


rmpath(CodeDir);

clear;
 
