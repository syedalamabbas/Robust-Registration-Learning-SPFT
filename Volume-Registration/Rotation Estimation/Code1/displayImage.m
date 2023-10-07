%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D:\Eye gaze\Spharm3DModel\3Ddataset\dataTest\test  XYZ positive -xy-z
%negative objectxyz  object-xy-z
clc

clear all;
close hidden all;
clear hidden all;
count=0;
    jjj=1;
InDataDir = '..\';
 DataFolder = 'benchmark\db\';
  DataFolder1 = 'benchmark\db\objectname\';
  
 for kkk=1:1
     
   for jj=1:100    
%inFiles = dir([InDataDir   DataFolder  num2str(kkk-1) '\m'  num2str(100*(kkk-1)+(jj-1)) '\m' num2str(100*(kkk-1)+(jj-1)) '_thumb.jpg']);
inFiles = dir([InDataDir   DataFolder1     '*.jpg']);

  inNames_obj={};
     for ii=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir   DataFolder1   inFiles(ii).name];
     end

    %*******************************************************************  
%% (1) Input data directory     
    for i=1:length(inFiles)
     %load(inNames_obj{i});
     
     I{jj}=imread(inNames_obj{i}); % Load the image file and store it as the variable I.


      [pa,name,ex]=fileparts(inNames_obj{i});     
      filetemplate=name;
   


 %PlotSurface1(vertices,faces);
    end
   end
       kkk
jjj=1;
   for m =1:5
       figure(m);
     for n=1:20
     subplot(4,5,n), imshow(imresize(I{jjj},4));
          title([num2str(kkk-1), 'class',num2str(jjj-1)]);
     jjj=jjj+1;
    %  hold on;

     end
          %title(['class',num2str(kkk-1),num2str(jjj)]);
    %close all;

   end


  

 end

