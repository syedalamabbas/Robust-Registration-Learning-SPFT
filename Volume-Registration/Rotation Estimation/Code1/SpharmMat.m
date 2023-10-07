

function [rmsd] = SpharmMat(confs, objs, method)

numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end
outNames = {};

%class(confs) - struct
%confs.vars -list of variable
%class(objs) - cell

h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file = objs{i};
    [path, name, ext, ver] = fileparts(file);
    
    switch method
        case 'Eigenvector'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_Eigenvector.log']));
            % Call AligSHREC method
            if isempty(confs.Template)
                disp('Eigenvector needs a template object');
                return;
            end
        
                [rmsd(1,i)] = eigenvector(file, confs);
                
               % [rmsd, outNames{end+1}] = align_CPS_SHREC(file, confs);
             %[vertices, sph_verts, faces, fvec, outNames{end+1}, confs] = align_CPS_SHREC(file, confs);
            clear('vertices', 'sph_verts', 'faces', 'fvec');
        
    end
    diary('off');
    waitbar(i/numSbj)
end
close(h);

return