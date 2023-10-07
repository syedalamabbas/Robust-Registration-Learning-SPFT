%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outNames = SpharmMatExpansion(confs, objs, method)

numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end
outNames = {};

h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file = objs{i};
    [path, name, ext] = fileparts(file);
    
    switch method
        case 'ExpLSF'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_ExpLSF.log']));
            % Call LSF method
            maxDeg = confs.MaxSPHARMDegree;
            [fvec, deg, Z, outNames{end+1}]=create_SPHARM_des_LSF([],[],[],maxDeg,file,deblank(char(confs.OutDirectory)));
              length(fvec)
            % fvecc=abs(fvec);
            %  save ('-ascii', ['..\3DdatasetFace\fvec.txt'], 'fvecc');
%               fid = fopen('..\3DdatasetFace\fvec.txt', 'w');
%              fprintf(fid, ' %16.2f %16.2f %16.2f \n', fvec);
%             fclose(fid);
            clear('fvec','fvecc','deg','Z');
    end
    diary('off');
    waitbar(i/numSbj)
    close all;
  %  clc;
end
close(h);

return