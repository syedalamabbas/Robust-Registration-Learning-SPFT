
function outNames = SpharmMatAlignment(confs, objs, method)

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
    [path, name, ext] = fileparts(file);
    
    switch method
        case 'AligSHREC'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_AlignSHREC.log']));
            % Call AligSHREC method
            if isempty(confs.Template)
                disp('SHREC needs a template object');
                return;
            end
           [vertices, sph_verts, faces, fvec, outNames{end+1}]=align_CPS_SHREC(file, confs);
            %clear('vertices', 'sph_verts', 'faces', 'fvec');
        case 'AligFOE'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_AlignFOE.log']));            
            % Call AligFOE method
            [vertices, sph_verts, faces, fvec, outNames{end+1}] = align_FOE(file, confs);
            clear('vertices', 'sph_verts', 'faces', 'fvec');
        case 'AligLandmark'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_AlignLDMK.log']));            
            % Call AligLandmark method
            disp('Not Implmented');
    end
    diary('off');
    waitbar(i/numSbj)
end
close(h);

return