

function outNames = SpharmMatParameterization(confs, objs, method)

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
        case 'ParamCALD'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_ParamCALD.log']));
            % Call CALD method (Intial parameterization and smoothing)
            [vertices_new, faces_new, sph_verts_new, outNames{end+1}] = SpharmMatParamTriaMesh(file, confs);
            clear('vertices_new', 'faces_new', 'sph_verts_new');
       
    end
    diary('off');
    waitbar(i/numSbj)
    close all;
    clc;
end
close(h);

return