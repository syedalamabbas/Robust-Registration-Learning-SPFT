clc
clear all;
close hidden all;
clear hidden all;
count=0;
InDataDir = '..\';
CodeDir = [InDataDir 'code'];
addpath(CodeDir);

%  DataFolder = 'benchmark\db\'; %object\test';
 DataFolder = '3Ddataset\object\test';
Templatenamedir=[InDataDir DataFolder '\Template\' ];

k=34;   
%data=zeros(115,5);
  n=75;   % n=25 hole 1.5 '(6+,18)'
  deg1=15;
  
  
 for kkk=6:6
   for jj=2:5 %40 19-24car
      kkk
     close hidden all; 
    
%     inFiles = dir([InDataDir   DataFolder  num2str(kkk-1) '\m'  num2str(100*(kkk-1)+(jj-1)) '\*.off']);
%     inNames_obj={};
%      for ii=1:length(inFiles)
%           inNames_obj{end+1} = [InDataDir   DataFolder  num2str(kkk-1) '\m'  num2str(100*(kkk-1)+(jj-1)) '\' inFiles(ii).name];
%      end

%      inNames_obj={};
%       ObjectName1=[InDataDir   DataFolder  num2str(kkk-1) '\Template\'  'ObjectName.txt'];
%     [ObjectName] = textread(ObjectName1, ' %s ' );   
%     for ii=1:length(ObjectName)
%          inNames_obj{end+1} = [InDataDir   DataFolder  num2str(kkk-1) '\'  ObjectName{ii}  '\' ObjectName{ii} '.off'];
%      end



  inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
     for ii=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
     end

%% (1) Input data directory     
   for i=1:length(inFiles)
 % for i=1:length(ObjectName)

 limt=4;  %34;
       load(inNames_obj{i});
 %     [vertices,faces] = read_mesh(inNames_obj{i});
      vertices=vertices';
      faces=faces';
      [pa,name,ex]=fileparts(inNames_obj{i});     

     origin=[0 0 0];
     vxsize =[1 1 1];
     confs.Connectivity='(6+,18)'; %
     confs.Epsilon = 1.5;
     confs.OutDirectory = [InDataDir DataFolder '/dataFix'];
     confs.name=name;
      con.name=name;
      con.flag=1;
     if ~exist(confs.OutDirectory,'dir')
         mkdir(confs.OutDirectory);
     end
    %///////////////////////////////////////////////////////////////////////////////////////////////////
    [bim] = verticestovolumefunc(vertices,faces,n);  
    
    [bim, new_name] = fix_holes(bim, confs);

      
%% (1) Input data directory     
   
    krow=1; % matrix to be eriten in the file
  %for k=0:5:limt  
     k=0;
      close all;
     % for R= 0:10:180
         R=0;
       [vertices, faces] =  gen_surf_data(bim,origin,vxsize); 
     % if (abs((length(vertices(:,1))-length(faces(:,1))) > 1))

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

%% (4.1) Perform Surface Meshes (CALD)

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

fvc33.vertices=vertices;
fvc33.faces=faces;
%[sph_verts,con,new_name] = initParamCALD1(fvc33,con);
 [sph_verts, new_name] = initParamCALD(vertices, faces, name, confs);


 %[vertices, faces, sph_verts, new_name] = smootheCALD(vertices, faces, sph_verts, new_name, confs);
%outNames = SpharmMatParameterization(confs, inNames, 'ParamCALD');
%///////////////////////////////inNames_obj   no rotation
if (con.flag==1)
%% (3) Display Surface Meshes (optional)

    %Available values for Space- 'object';'param';'both'
dispConfs.Space = 'both';
    % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    %    'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'orig';  
    % Available values for Shape- 'solid';'mesh';'both'
dispConfs.Shade = 'both';
    % Available values for Overlay- 'none';'adc_paramap'
dispConfs.Overlay = 'adc_paramap';
    % Available values for Export- 'screen';'png';'both'
dispConfs.Export = 'both';
dispConfs.Degree = [];
dispConfs.Template = '';
dispConfs.count=R;
dispConfs.name=con.name;
SpharmMatUtilDisplayObjsParam(dispConfs, new_name, CodeDir);

    %  end %R rotation
 %  end %limit;
   close all;
end
 end % length of file
 
    % end %if
%   if((sum==0) && (L >0))
%       Objectname= [ Templatenamedir   'ObjectName.txt'];
%      fid = fopen(Objectname, 'a');
%      fprintf(fid, ' %s  \n', name );
%      fclose(fid);
%   end
        
      

     end  % length of the objects  jj
 end %length of the class  kkk
