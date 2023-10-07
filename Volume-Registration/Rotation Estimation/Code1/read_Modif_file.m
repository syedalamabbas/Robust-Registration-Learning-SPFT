clc

clear all;
close hidden all;
clear hidden all;
ii=0;
%function [data, zmin, nrows, ncols, imfile] = read_bntfile(filepath)
for i=0:104
InDataDir = '..\'; 

if((i<10) && (i>=0))
 DataFolder = ['\RResult75\' 'bs00' num2str(i)  ]; 
 inFiles=dir([InDataDir   DataFolder '\bs00' num2str(i) '_E_*' '*.mat']);
elseif ((i>9) && (i<100))
  DataFolder = ['\RResult75\' 'bs0'  num2str(i)  ]; 
 inFiles=dir([InDataDir   DataFolder '\bs0' num2str(i) '_E_*' '*.mat']);   
elseif (i>99) 
 DataFolder = ['\RResult75\' 'bs'  num2str(i)  ]; 
 inFiles=dir([ InDataDir  DataFolder '\bs' num2str(i) '_E_*' '*.mat']);   
end  %end if





if (length(inFiles)==6)
  inNames_obj={};
     for j=1:length(inFiles)
          inNames_obj{end+1} = [InDataDir  DataFolder '\' inFiles(j).name];
     end
     
for j=1:length(inFiles)
    load( inNames_obj{j});
 
if((ii<10) && (ii>=0))
  [pa,name,ex,ve]=fileparts(inNames_obj{j});
 DataFoldersave = ['\Result75\' 'bs00' num2str(ii)]; 
 OutDirectory = [InDataDir DataFoldersave ];
  if ~exist(OutDirectory,'dir')
      mkdir(OutDirectory);
  end
   DataFoldersave = ['\Result75\' 'bs00' num2str(ii) '\bs00' num2str(ii) name(6:end) '.mat']; 
   OutDirectory = [InDataDir DataFoldersave ];


elseif ((ii>9) && (ii<100))
 [pa,name,ex,ve]=fileparts(inNames_obj{j});
 DataFoldersave = ['\Result75\' 'bs0'  num2str(ii)]; 
 OutDirectory = [InDataDir DataFoldersave ];
  if ~exist(OutDirectory,'dir')
      mkdir(OutDirectory);
  end
     DataFoldersave = ['\Result75\' 'bs0' num2str(ii) '\bs0' num2str(ii) name(6:end) '.mat' ]; 
   OutDirectory = [InDataDir DataFoldersave ];
   
end  %end if
    save(OutDirectory, 'vertices', 'faces', 'sph_verts', 'fvec');
end  %end for
ii=ii+1;
end % end if length(inFiles)

end %for i
