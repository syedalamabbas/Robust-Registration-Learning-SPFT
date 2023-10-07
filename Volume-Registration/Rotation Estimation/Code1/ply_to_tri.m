%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D:\Eye gaze\Spharm3DModel\3Ddataset\dataTest\test  XYZ positive -xy-z
%negative objectxyz  object-xy-z
clc

clear all;
close hidden all;
clear hidden all;

InDataDir = '..\';
 DataFolder = 'benchmark\missingOBJECT\test';
 DataFolder1 = 'benchmark\missingOBJECT\test\';
  inFiles = dir([InDataDir   DataFolder '\*.ply']); 
  inNames_obj={};
     for i=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir   DataFolder '\' inFiles(i).name];
     end

 for i=1:length(inFiles)
  %% (1) Input data directory     
      [vertices,faces] = read_mesh(inNames_obj{i});
      [pa,name,ex]=fileparts(inNames_obj{i});
      
     % fileN=name(1:end-4);
       fileN=name;
      ply_file_name=[InDataDir   DataFolder1 fileN '.ply'];
        confs.OutDirectory = [InDataDir DataFolder1];     
      if ~exist(confs.OutDirectory,'dir')
             mkdir(confs.OutDirectory);
      end
      % plotsurfaceT(vertices,faces,filetemplate,theta,theta_r);
      % PlotSurface1(vertices,faces);
        tri_file_name=[InDataDir   DataFolder1 fileN '.mat'];
      
        save(tri_file_name,'vertices','faces');

          i
 end

