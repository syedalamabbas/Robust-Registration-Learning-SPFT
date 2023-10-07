

function [rmsd, theta] = SpharmMethodtheta(confs, objs, method, theta)

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
            % Call Eigenvector method
            if isempty(confs.Template)
                disp('Eigenvector needs a template object');
                return;
            end
               [rmsd(1,i)] = eigenvector(file, confs);
                 break;
        case 'FOE'
            % Call FOE method
             [rmsd(1,i)] = FOE(file, confs);
                 break;
                 
       case 'ICP_Eigenvector'
            % Call Eigenvector method
            if isempty(confs.Template)
                disp('ICP_Eigenvector needs a template object');
                return;
            end
               [rmsd(1,i)] = ICP_eigenvector(file, confs);
                 break;
       case 'ICP'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_AlignLDMK.log']));            
            % Call ICP method
            disp('Not Implmented');
            break;
            
       otherwise
              disp('No method defined.');
              rmsd=0;
    end  
    
      diary('off');
    waitbar(i/numSbj)
end
close(h);

return