

function outNames = SpharmMatAlignment1(confs, objs)

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
    [path, name, ext, ver] = fileparts(file);
    
         % Call AligFOE method
            [vertices, sph_verts, faces, fvec, outNames{end+1}] = align_FOE1(file, confs);
            clear('vertices', 'sph_verts', 'faces', 'fvec');
       
    waitbar(i/numSbj)
end
close(h);

return