function Ne=vertex_neighbours(vertices, faces)
% This function VERTEX_NEIGHBOURS will search in a face list for all 
% the neigbours of each vertex.
%
% Ne=vertex_neighbours(FV)
%

Ne=vertex_neighbours_double(faces(:,1),faces(:,2),faces(:,3),vertices(:,1),vertices(:,2),vertices(:,3));
