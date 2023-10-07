clc
clear all;
close hidden all;
clear hidden all;
count=0;
InDataDir = '..\';

DataFolder = 'benchmark\3DObject_Dataset\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];
k=7;
deg1=15;


close hidden all;


inFiles = dir([InDataDir   DataFolder '\*.mat']); inNames_obj={};
for ii=1:length(inFiles)
    inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(ii).name];
end

fid = fopen('Result.csv', 'a');
fprintf(fid, '\n Object Name ,  RMSD , RMSD0  , Time(sec), Delta_x, Delta_y, Delta_z, Actual_x , Actual_y , Actual_z ,Norm_Difference \n');
fclose(fid) ;


%% (1) Input data directory
for i=1:9; %length(inFiles)
    %  for i=1:length(ObjectName)
    
    limt=4;  %35;
    load(inNames_obj{i});
    [pa,name,ex]=fileparts(inNames_obj{i});
    
    
    PlotSurface1(vertices,faces);
    
    
    origin=[0 0 0];
    vxsize =[1 1 1];
    confs.Connectivity='(6+,18)'; 
    confs.Epsilon = 5.5000;
    confs.OutDirectory = [InDataDir DataFolder 'dataFix'];
    confs.name=name;
    if ~exist(confs.OutDirectory,'dir')
        mkdir(confs.OutDirectory);
    end
    %///////////////////////////////////////////////////////////////////////////////////////////////////
    [bim] = verticestovolumefunc(vertices,faces);
     
    [bim, new_name] = fix_holes(bim, confs);
    
    
    %% (1) Input data directory
    
    krow=1; % matrix to be eriten in the file
    for k=0:5:limt
        close all;
        for R=[23,34.56,57.78,87.34]
            
%             N = 55;
%             maxVal = 24;
%             if(R ~= 0)  
%                 for s = 1:N
%                     J = bim(:,:,s);
%                     randomIntegerL = randi(maxVal);
%                     randomIntegerR = randi(maxVal);
%                     J((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR,(N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR) = zeros(length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR),length((N-1)/2-randomIntegerL :(N-1)/2+ randomIntegerR));
%                     bim(:,:,s) = J;
%                 end
%                 % Before Hole Filling
%                  [vertices, faces] =  gen_surf_data(bim,origin,vxsize);
%                 PlotSurface1(vertices,faces);
%                 
%                 [bim, new_name] = fix_holes(bim, confs);   % It fails with or without this
%             end
                    
            
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
            if(R == 0)
                theta_r.x=R*kx;
                theta_r.y=R*ky;
                theta_r.z=R*kz;
            else
                theta_r.x=R*kx;
                theta_r.y=(34+rand(1,1))*ky;
                theta_r.z=(64+rand(1,1))*kz;
            end
%             disp(sprintf('Ture rotation  [theta]xyz: %0.2f %0.2f %0.2f',(theta_r.x),(theta_r.y),(theta_r.z)));
            
            
            [vertices,blob_center] = center(vertices);
            [vertices] = rotate_obj3d(vertices,theta_r);
            PlotSurface1(vertices,faces);
            
            tic
             
            %% (4.1) Perform Surface Meshes (CALD)
            
            %Create spherical parameterization for surface meshes  **Mapping**
            
%             deg1 = 31;
            
            confs.MeshGridSize = 50;
            confs.MaxSPHARMDegree = deg1;
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
            
            %% (3) Display Surface Meshes (optional)
            
                %Available values for Space- 'object';'param';'both'
%             dispConfs.Space = 'both';
%                 % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
%                 %    'icosa3';'icosa4';'icosa5';'icosa6'
%             dispConfs.Mesh = 'orig';
%                 % Available values for Shape- 'solid';'mesh';'both'
%             dispConfs.Shade = 'both';
%                 % Available values for Overlay- 'none';'adc_paramap'
%             dispConfs.Overlay = 'adc_paramap';
%                 % Available values for Export- 'screen';'png';'both'
%             dispConfs.Export = 'both';
%             dispConfs.Degree = [];
%             dispConfs.Template = '';
%             dispConfs.count=R;
%             SpharmMatUtilDisplayObjs(dispConfs, new_name, CodeDir);
            %
            
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
            
            Templatename=[Templatenamedir   name  '_tmp.mat'];
            if (R==0)
                if ~exist(Templatenamedir,'dir')
                    mkdir(Templatenamedir);
                end
                save(Templatename,'vertices','faces','fvec');
            end
           
            
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
            
            
            
            vertices1=vertices;
            load(confs.Template);
            vertices2=vertices;
            vertices2=vertices2*RR;
            rmsdO=((vertices1-vertices2));
            rmsdO=norm(rmsdO(:))/sqrt(length(vertices2(:,1)));
            %rmsd1=(norm (abs(vertices1-vertices2)))
            
            t= toc;
            disp(sprintf('Estimated rotation  [theta]xyz: %0.2f %0.2f %0.2f',(theta.x),(theta.y),(theta.z)));
            
            data(krow,1)=rmsd;
            data(krow,2)=rmsdO;
            data(krow,3)=t;
            data(krow,4)=theta.x-theta_r.x;
            data(krow,5)=theta.y-theta_r.y;  
            data(krow,6)=theta.z-theta_r.z;
            data(krow,7)=theta_r.x;
            data(krow,8)=theta_r.y;
            data(krow,9)=theta_r.z;
            
            
            %% Write to file 
            
            normVal = norm(cell2mat(struct2cell(theta_r)) - cell2mat(struct2cell(theta)));
            fid = fopen('Result.csv', 'a');
            fprintf(fid, '\n %s ,  %2.2f , %2.2f  , %2.2f ,%2.2f, %2.2f, %2.2f, %2.2f , %2.2f ,  %2.2f , %2.9f ',  name, data(krow,:), normVal );
            fclose(fid) ;
            
            krow=krow+1;
        end %R rotation
        %               PlotSurface1(vertices,faces);
    end %limit;
    
    close all;
    clear data
end % length of file
