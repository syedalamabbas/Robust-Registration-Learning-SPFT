

function [vertices, faces] = readM(filename)

run(filename);

% If p is the handle of patch, created from the above command
% run(fileListOnDisk(1,:));
vertices = get(p, 'Vertices');
faces = get(p, 'Faces');

return;