
function  [rmsd] = SpharmMatAlig(confs, objs)

numSbj = length(objs);

for i = confs.count:numSbj
    file = objs{i};
         [rmsd(1,i)] = alig(file, confs);  

end

return