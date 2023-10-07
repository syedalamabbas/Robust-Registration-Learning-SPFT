

function [names] = Compare(confs, objs)

numSbj = length(objs);


h = waitbar(0,'Please wait...');
for i = confs.count: numSbj
    file = objs{i};
    [path, name1, ext, ver] = fileparts(file);
        [path, name2, ext, ver] = fileparts(confs.Template);
    
 
           % diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_Eigenvector.log']));
            % Call String Comapre
            
            R=strncmpi(name1, name2', 9);
            
            if isempty(confs.Template)
                disp('Eigenvector needs a template object');
                return;
            end

            if (R==1)
                names{i}= name1;
                        
            end
            
           
   
    diary('off');
    waitbar(i/numSbj)
end
close(h);

return