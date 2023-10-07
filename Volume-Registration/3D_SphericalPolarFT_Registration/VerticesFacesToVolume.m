function [ Volume ] = VerticesFacesToVolume( vertices,faces, N )
%VERTICESFACESTOVOLUME Summary of this function goes here
% change quadralaterals to triangles

% faces = delaunay(vertices(:,1), vertices(:,2), vertices(:,3));

if size(faces,2)==4
    qfaces = faces;
    faces = [faces(:,1:3); faces(:,[3 4 1])];
    
    dif1 = faces(:,1)-faces(:,2); dif2=faces(:,2)-faces(:,3); dif3=faces(:,3)-faces(:,1);
    indDif1 = find(dif1 == 0);
    indDif2 = find(dif2 == 0);
    indDif3 = find(dif3 == 0);
    indDif = union([indDif1, indDif2, indDif3],'legacy');
    ufIDX = setdiff([1:size(faces,1)], indDif);
    faces = faces(ufIDX,:);
end% data=vertices .*50;

fvc3.vertices=vertices;
fvc3.faces=faces;
% N=55;
  % Convert the mesh to a voxelvolume
  Volume=polygon2voxel(fvc3,[N N  N],'auto');
 % Volume(:,:,1:10)=0;

end

