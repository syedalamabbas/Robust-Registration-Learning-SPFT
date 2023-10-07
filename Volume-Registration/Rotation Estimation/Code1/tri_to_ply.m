%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D:\Eye gaze\Spharm3DModel\3Ddataset\dataTest\test  XYZ positive -xy-z
%negative objectxyz  object-xy-z
clc

clear all;
close hidden all;
clear hidden all;

InDataDir = '..\';
 DataFolder = 'benchmark\db\objectname\';
 DataFolder1 = 'benchmark\db\objectname\';
  inFiles = dir([InDataDir   DataFolder '\*.mat']); 
  inNames_obj={};
     for i=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(i).name];
     end

 for i=383:383 %length(inFiles)
  %% (1) Input data directory     
      load(inNames_obj{i});
      [pa,name,ex]=fileparts(inNames_obj{i});
      
     % fileN=name(1:end-4);
       fileN=name;
      filename=[InDataDir   DataFolder1 fileN '.ply'];
        confs.OutDirectory = [InDataDir DataFolder1];     
      if ~exist(confs.OutDirectory,'dir')
             mkdir(confs.OutDirectory);
      end
      % plotsurfaceT(vertices,faces,filetemplate,theta,theta_r);
      % PlotSurface1(vertices,faces);
      
qfaces=[];

% change quadralaterals to triangles
if size(faces,2)==4
    qfaces = faces;
    faces = [faces(:,1:3); faces(:,[3 4 1])];
    
    dif1 = faces(:,1)-faces(:,2); dif2=faces(:,2)-faces(:,3); dif3=faces(:,3)-faces(:,1);
    indDif1 = find(dif1 == 0);
    indDif2 = find(dif2 == 0);
    indDif3 = find(dif3 == 0);
    indDif = union(indDif1, indDif2, indDif3);
    ufIDX = setdiff([1:size(faces,1)], indDif);
    faces = faces(ufIDX,:);
end

      
      
      
       ply_data = tri_surface_to_ply ( vertices', faces' );
       
       ply_write(ply_data,filename,'ascii');
       i


 end

