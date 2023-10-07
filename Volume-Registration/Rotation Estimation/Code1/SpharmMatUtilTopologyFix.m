
function outNames = SpharmMatUtilTopologyFix(confs, objs, method)
%numSbj = size(objs,2);
numSbj = 1;
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

outNames = {};

h = waitbar(0,'Please wait...');
for i = 1:numSbj
   file = objs;
   % file = objs{i};    
    [pa,na,ex]=fileparts(file);
    
    switch method
        case 'InHouse_Fix'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' na '_fix_InHouse.log']));
            [bim_new,origin_new,vxsize_new,outNames{end+1}] = fix_bad_topology(file, confs);
        
           
    end
    diary('off');
    waitbar(i/numSbj)
    close all;
end
close(h);

return;
