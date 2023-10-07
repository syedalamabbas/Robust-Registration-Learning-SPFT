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

 DataFolder = 'benchmark\3DObject_Dataset\objectname\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];

  deg1=15;
 class={ 'Flower'  'Car'};
     limt=6;  %35;
 for c=1:2
       Templatenamedir=[InDataDir DataFolder   class{c}  '\dataOutput\' ];
  for object=1:5 

  %inFiles = dir([InDataDir   DataFolder   class{c}  '\dataOutput\'   '*FOE_als.mat']);
  inFiles = dir([InDataDir   DataFolder   class{c}  '\dataExp\'   '*.mat']);
  inNames_obj={};
     for ii=1:length(inFiles);
          
        inNames_obj{end+1} = [InDataDir   DataFolder  class{c}  '\dataExp\'  inFiles(ii).name];
     end

    %*******************************************************************  
%% (1) Input data directory     
     for i=1:length(inFiles);  
     load(inNames_obj{i});
%       [vertices,faces] = read_mesh(inNames_obj{i});
%       vertices=vertices';
%       faces=faces';
      [pa,name,ex]=fileparts(inNames_obj{i});     
      filetemplate=name;
      origin=[0 0 0];
      vxsize =[1 1 1];
      confs.Connectivity='(6+,18)'; 
      confs.Epsilon = 1.5000;
      confs.OutDirectory = [InDataDir DataFolder '/dataFix'];
      confs.name=name;
       if ~exist(confs.OutDirectory,'dir')
              mkdir(confs.OutDirectory);
       end

        [bim] = verticestovolumefunc(vertices,faces,n);

        [bim, new_name] = fix_holes(bim, confs);

        [vertices, faces] =  gen_surf_data(bim,origin,vxsize);
        
         sum=0;

    krow=1; % matrix to be eriten in the file
      for k=0:5:limt        
         for R= 0:10:80
         
   
           close all;
          if (k==0)
               kz=0;
               ky=0;
               kx=0;  
          elseif (k==5)
             % around x y Z
              kz=1;
              ky=1;
              kx=1;
         elseif (k==10)
               %       around x
             kz=0;
             ky=0;
             kx=(k-9);
         elseif(k==15)
          %    around Z
            kz=(k-14);
            ky=0;
            kx=0;
          elseif (k==20)
          %    around Y
           kz=0;
           ky=(k-19);
           kx=0; 
    
         elseif(k==25)
           %      around X Y
           kz=0;
           ky=(k-24);
           kx=(k-24);

        elseif (k==30)
          % around X Z
         kz= (k-29);
         ky= 0;
         kx= (k-29);

        elseif (k==35)
        %  around Y  Z
        kz=(k-34);
        ky=(k-34);
        kx=0;
        else
        kz=0;
        ky=0;
        kx=0;
        end


    [vertices, faces] =  gen_surf_data(bim,origin,vxsize);
%      vertnum = size(vertices,1);
%      maxDeg= max(1, floor(sqrt(vertnum)*1/2));
     maxDeg=deg1;
    %% % (4) Make 3D Rotation
    theta_r.x=R*kx;
    theta_r.y=R*ky;
    theta_r.z=R*kz;
     [vertices,blob_center] = center(vertices);
    [vertices] = rotate_obj3d(vertices,theta_r);
    
%% (4.1) Perform Surface Meshes (CALD)

%Create spherical parameterization for surface meshes  **Mapping**
%// Read Oject
confs.MeshGridSize = 50;
confs.MaxSPHARMDegree = maxDeg;  % 6
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

tic;
[sph_verts, new_name] = initParamCALD(vertices, faces, name, confs);

 %[vertices, faces, sph_verts, new_name] = smootheCALD(vertices, faces, sph_verts, new_name, confs);
%outNames = SpharmMatParameterization(confs, inNames, 'ParamCALD');
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
% 
% SpharmMatUtilDisplayObjs(dispConfs, new_name, CodeDir);
% 

%% (5) Perform SPHARM-MAT Expansion
% Calculate SPHARM coefficients
confs.MaxSPHARMDegree = maxDeg;   % 15
confs.OutDirectory = [InDataDir DataFolder '/dataOutput' '/dataExp'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
outNames = {};
 %outNames = SpharmMatExpansion(confs, outNames, 'ExpLSF');
 [fvec, deg, Z, outNames{end+1}]=create_SPHARM_des_LSF([],[],[],maxDeg,new_name,deblank(char(confs.OutDirectory)));
 
%   Templatename=[Templatenamedir    'atlas.mat'];
% if (R==0)
%   if ~exist(Templatenamedir,'dir')
%     mkdir(Templatenamedir);
%   end  
%   save(Templatename,'vertices','faces','fvec');
% end

%% (7)estimate Rotation

confs.Template = [Templatenamedir   'atlas.mat'];  % obj file same file

confs.MaxSPHARMDegree = maxDeg;   %12   depends on the vertices number deg = max(1, floor(sqrt(vertnum)*1/2));
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
% if (theta.x >85) 
%     theta.x=thetaxp;
% end
% if (theta.y >85) 
%     theta.y=thetayp;
% end
% if (theta.z >85) 
%     theta.z=thetazp;
% end

vertices1=vertices;
load(confs.Template);
 vertices2=vertices;
 vertices2=vertices2*RR;
rmsdO=((vertices1-vertices2));
rmsdO=norm(rmsdO(:))/sqrt(length(vertices2(:,1)))
%rmsd1=(norm (abs(vertices1-vertices2)))

t= toc;
  disp(sprintf('Estimated rotation  [theta]xyz: %0.2f %0.2f %0.2f',(theta.x),(theta.y),(theta.z)));
% sum =sum+(theta.x-R)+(theta.y-R)+(theta.z-R);
%   if(sum > 10)
%      break;
%   end
data(krow,1)=rmsd;
data(krow,2)=rmsdO;
data(krow,3)=t;
data(krow,4)=theta.x-theta_r.x;
data(krow,5)=theta.y-theta_r.y;
data(krow,6)=theta.z-theta_r.z;
data(krow,7)=theta_r.x;
data(krow,8)=theta_r.y;
data(krow,9)=theta_r.z;
thetaxp=theta.x;
thetayp=theta.y;
thetazp=theta.z;
%[rmsd] = SpharmMatAlignment(confs, outNames, 'AligSHREC');
%[ rmsd, confs] = SpharmMat(confs, outNames, 'Eigenvector');   


       
%% Write to file
%data= data(find(data(:,1)~= 0),:);

% Rx = corr2(datax(:,1),datax(:,2))
% Ry = corr2(datay(:,1),datay(:,2))
% Rz = corr2(dataz(:,1),dataz(:,2))
% meanx=mean(data(:,2));  % Mean of the abs error
% %  medianx=median(datax(:,3));
% stdx=std(data(:,2));   % STD of the absoulte error
% % sumx=sum(datax(:,3));   % Sum of the absoulte error
% 
% meany=mean(data(:,3));
% %  mediany=median(datay(:,3));
% stdy=std(data(:,3));
% % sumy=sum(datay(:,3));
% 
% meanz=mean(data(:,4));
% %  medianz=median(dataz(:,3));
% stdz=std(data(:,4));  
% % sumz=sum(data(:,3));
% M_rmsd=mean(data(:,1));  % RMSD
% % M_minutes=mean(data(:,5));   % seconds;
% % M_rmsdO=mean(data(:,6));  % RMSDO
% % [VROW,VCOL]= size(vertices);
% 
% 
% %    filetxtname=['..\3Ddataset\datatest\test1\pngR\' name  '.txt'];
% %      save ('-ascii', ['..\3Ddataset\datatest\test\pngR\' name  'x.txt'], 'datax');
% %      save ('-ascii', ['..\3Ddataset\datatest\test\pngR\' name  'y.txt'], 'datay');
% %      save ('-ascii', ['..\3Ddataset\datatest\test\pngR\' name  'z.txt'], 'dataz');
% %      save ('-ascii', ['..\3Ddataset\datatest\test\pngR\' name  'TL.txt'], 'dataTL');
% %      
% %      fid = fopen('..\3Ddataset\datatest\test\Result_XYZ\Mean_Data.txt', 'a');
%       fid = fopen('..\benchmark\db\Result\Mean_Data.txt', 'a');
%      mean_xyz=[meanx, meany,meanz];
%      mean_xyz= mean(mean_xyz);
%      fprintf(fid, ' %s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', name, meanx, meany,meanz,mean_xyz, theta_r.x, theta_r.y, theta_r.z);
%      fclose(fid);
% %      fid = fopen('..\3Ddataset\datatest\test\Result_XYZ\STD_Data.txt', 'a');
%      fid = fopen('..\benchmark\db\Result\STD_Data.txt', 'a');
%      STD_xyz=[stdx,stdy,stdz];
%      STD_xyz= std(STD_xyz);
%      fprintf(fid, ' %s %6.2f %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f\n', name, stdx,stdy,stdz,STD_xyz, theta_r.x, theta_r.y, theta_r.z);
%      fclose(fid);

  
     
%       fid = fopen('..\benchmark\db\Result\Data.txt', 'a');
%      fprintf(fid, '%s  %4.2f %4.2f  %4.2f   %4.2f %4.2f %4.2f %4.2f %4.2f  %4.2f\n',  name, data(krow,:) );

%      fclose(fid) ;
 datarmsdO(1,krow)=rmsdO;
 datarmsd(1,krow)=rmsd;
     krow=krow+1;
          end %R rotation
       end %limit;
   
      fid = fopen('..\benchmark\db\Result\DatarmsdO.txt', 'a');
     fprintf(fid, '\n %s  %5f %2.4f  %2.4f   %2.4f %2.4f %2.4f %2.4f %2.4f  %2.4f\n',  name,length(v(:,1)), datarmsdO(krowR,:) );

     fclose(fid) ;  
     fid = fopen('..\benchmark\db\Result\DatarmsdF.txt', 'a');
     fprintf(fid, '\n %s  %5f %2.4f  %2.4f   %2.4f %2.4f %2.4f %2.4f %2.4f  %2.4f\n',  name,length(v(:,1)), datarmsd(krowR,:) );

     fclose(fid) ;
       end % length of file
 
 %    end %if
%   if((sum==0) && (L >0))
%       Objectname= [ Templatenamedir   'ObjectName.txt'];
%      fid = fopen(Objectname, 'a');
%      fprintf(fid, ' %s  \n', name );
%      fclose(fid);
%   end
        
      

  end  % length of the objects
 end %length of the class
