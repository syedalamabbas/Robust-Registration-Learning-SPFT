clc
clear all;
close hidden all;
clear hidden all;
count=0;
InDataDir = '..\';
CodeDir = [InDataDir 'code'];
addpath(CodeDir);

 DataFolder = 'benchmark\db\';
Templatenamedir=[InDataDir DataFolder '\Template\' ];

k=34;   
%data=zeros(115,5);
  n=25;
  deg1=15;
  
  
 for kkk=17:17  
  %for jj=1:99 %40 19-24car
      kkk
     close hidden all; 
    
%     inFiles = dir([InDataDir   DataFolder  num2str(kkk-1) '\m'  num2str(100*(kkk-1)+(jj-1)) '\*.off']);
%     inNames_obj={};
%      for ii=1:length(inFiles)
%           inNames_obj{end+1} = [InDataDir   DataFolder  num2str(kkk-1) '\m'  num2str(100*(kkk-1)+(jj-1)) '\' inFiles(ii).name];
%      end

     inNames_obj={};
      ObjectName1=[InDataDir   DataFolder  num2str(kkk-1) '\Template\'  'ObjectName.txt'];
    [ObjectName] = textread(ObjectName1, ' %s ' );   
    for ii=1:length(ObjectName)
         inNames_obj{end+1} = [InDataDir   DataFolder  num2str(kkk-1) '\'  ObjectName{ii}  '\' ObjectName{ii} '.off'];
     end
%% (1) Input data directory     
 %for i=1:length(inFiles) 
   for i=1:length(ObjectName)


           %load(inNames_obj{i});
      [vertices,faces] = read_mesh(inNames_obj{i});
      vertices=vertices';
      faces=faces';
      [pa,name,ex]=fileparts(inNames_obj{i});     

     origin=[0 0 0];
     vxsize =[1 1 1];
     confs.Connectivity='(6+,18)'; 
     confs.Epsilon = 5.5000;
     confs.OutDirectory = [InDataDir DataFolder '/dataFix'];
     confs.name=name;
     if ~exist(confs.OutDirectory,'dir')
         mkdir(confs.OutDirectory);
     end
    %///////////////////////////////////////////////////////////////////////////////////////////////////
    [bim] = verticestovolumefunc(vertices,faces,n);  
    
    [bim, new_name] = fix_holes(bim, confs);

      PlotSurface1(vertices,faces);
       
    [vertices, faces] =  gen_surf_data(bim,origin,vxsize); 
     
  end
 %end %jj
 end