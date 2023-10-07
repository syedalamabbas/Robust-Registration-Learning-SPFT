clc
clear all;
close hidden all;
clear hidden all;
count=0;
InDataDir = '..\';



 DataFolder = 'benchmark\3DObject_Dataset\3D_Object_Noise\'; 
Templatenamedir=[InDataDir DataFolder '\TemplateNoise\' ];

k=34;    
  deg1=15;
   
     inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
     for ii=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
     end
     

%% (1) Input data directory     
 for i=1:length(inFiles)
%  for i=1:length(ObjectName)

 limt=4;  %35;
           load(inNames_obj{i});
                 [pa,name,ex]=fileparts(inNames_obj{i}); 
                 
                 PlotSurface1(vertices,faces);

     origin=[0 0 0];
     vxsize =[1 1 1];
     confs.Connectivity='(6+,18)'; %
     confs.Epsilon = 1.5;
     confs.OutDirectory = [InDataDir DataFolder '/dataFix'];
     confs.name=name;
     if ~exist(confs.OutDirectory,'dir')
         mkdir(confs.OutDirectory);
     end
    %///////////////////////////////////////////////////////////////////////////////////////////////////
    %Adding noise //////////////////////////////////////////////////////////////////////////
     noise = randn(size(vertices))*0.1;
     vertices = vertices+  noise;   

%       
    [bim] = verticestovolumefunc(vertices,faces);  
    [bim, new_name] = fix_holes(bim, confs);


%% (1) Input data directory     
   
    krow=1; % matrix to be eriten in the file
  for k=0:5:limt  
      close all;
      for R= 0:10:80
     [vertices, faces] =  gen_surf_data(bim,origin,vxsize); 

         PlotSurface1(vertices,faces);  


        if (k==0)
             % around x y Z
              kz=1;
              ky=1;
              kx=1;
         elseif (k==5)
               %       around x
             kz=0;
             ky=0;
             kx=(k-4);
         elseif(k==10)
          %    around Z
            kz=(k-9);
            ky=0;
            kx=0;
          elseif (k==15)
          %    around Y
           kz=0;
           ky=(k-14);
           kx=0; 
    
         elseif(k==20)
           %      around X Y
           kz=0;
           ky=(k-19);
           kx=(k-19);

        elseif (k==25)
          % around X Z
         kz= (k-24);
         ky= 0;
         kx= (k-24);

        elseif (k==30)
        %  around Y  Z
        kz=(k-29);
        ky=(k-29);
        kx=0;
        else
        kz=0;
        ky=0;
        kx=0;
        end

    %% % (4) Make 3D Rotation
    theta_r.x=R*kx;
    theta_r.y=R*ky;
    theta_r.z=R*kz;
   %disp(sprintf('Ture rotation  [theta]xyz: %0.2f %0.2f %0.2f',(theta_r.x),(theta_r.y),(theta_r.z)));   
    
%     theta_r.x=40;
%     theta_r.y=80;
%     theta_r.z=80;


%      [fvec,blob_center] = center(fvec);
%     [fvec] = rotate_obj3d(fvec,theta_r);  % works with eigenvector not
    
    [vertices,blob_center] = center(vertices);
    [vertices,TR] = rotate_obj3d1(vertices,theta_r);
     V_True = vrrotmat2vec(TR);
tic

%% (4.1) Perform Surface Meshes 

%Create spherical parameterization for surface meshes  **Mapping**
%// Read Oject
confs.MeshGridSize = 50;
confs.MaxSPHARMDegree = deg1;  % 6
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


[sph_verts, new_name] = initParamCALD(vertices, faces, name, confs);

%///////////////////////////////inNames_obj   no rotation
%% (3) Display Surface Meshes (optional)

%     %Available values for Space- 'object';'param';'both'
% dispConfs.Space = 'both';
%     % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
%     %    'icosa3';'icosa4';'icosa5';'icosa6'
% dispConfs.Mesh = 'orig';  
%     % Available values for Shape- 'solid';'mesh';'both'
% dispConfs.Shade = 'both';
%     % Available values for Overlay- 'none';'adc_paramap'
% dispConfs.Overlay = 'adc_paramap';
%     % Available values for Export- 'screen';'png';'both'
% dispConfs.Export = 'both';
% dispConfs.Degree = [];
% dispConfs.Template = '';
% dispConfs.count=R;
% SpharmMatUtilDisplayObjs(dispConfs, new_name, CodeDir);
% % 

%% (5) Perform SPHARM-MAT Expansion
% Calculate SPHARM coefficients
confs.MaxSPHARMDegree = deg1;   % 15
confs.OutDirectory = [InDataDir DataFolder '/dataOutput' '/dataExp'];
maxDeg=deg1;
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
outNames = {};
 [fvec, deg, Z, outNames{end+1}]=create_SPHARM_des_LSF([],[],[],maxDeg,new_name,deblank(char(confs.OutDirectory)));


%% (7)estimate Rotation

confs.Template = [Templatenamedir  name  '_tmp.mat'];  % obj file same file
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
[ rmsd, theta,RR] = SpharmMethod(confs, outNames, method);  

V_Exp = vrrotmat2vec(real(RR));
V_True_unit=V_True(1,1:3)/norm(V_True(1,1:3));
V_Exp_unit=V_Exp(1,1:3)/norm(V_Exp(1,1:3));
% tt=rt(1,4)*180/pi;

V_Error= norm(V_True(1,1:3)-V_Exp(1,1:3));
A_Error=abs(V_True(1,4)-V_Exp(1,4));
A1_Error=acos(dot(V_True_unit,V_Exp_unit));

data(krow,1)=V_Error;
data(krow,2)=A_Error;
data(krow,3)=A1_Error;
data(krow,4)=theta.x;
data(krow,5)=theta.y;
data(krow,6)=theta.z;
data(krow,7)=theta_r.x;
data(krow,8)=theta_r.y;
data(krow,9)=theta_r.z;

  disp(sprintf('Estimated rotation  [theta]xyz: %0.2f %0.2f %0.2f',(theta.x),(theta.y),(theta.z)));
      
%% Write to file

% 
%       fid = fopen('..\benchmark\db1\Result\Resultdb13.txt', 'a');
%       fprintf(fid, '\n %s   %2.2f  %2.2f   %2.2f %2.2f %2.2f %2.2f %2.2f  %2.2f   %2.2f \n',  name, data(krow,:) );
%      fclose(fid) ;  
     
      krow=krow+1;
     % end  %if
      end %R rotation
   end %limit;
   close all;
 end % length of file
 

