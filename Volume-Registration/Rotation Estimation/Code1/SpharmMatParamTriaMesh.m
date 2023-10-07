%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vertices, faces, sph_verts, new_name] = SpharmMatParamTriaMesh(file, confs)

load(file);
[path, name, ext] = fileparts(file);

% First, if obj == bim, then extract surface
if exist('bim','var')
    % Make binary image
    ind = find(bim>0); bim(ind) = 1;
    ind = find(bim<1); bim(ind) = 0;    
    
    [vertices, faces] =  gen_surf_data(bim,origin,vxsize);
    %[faces,vertices]=voxel_bnd_faces(bim,vxsize,origin,1);
    
    postfix = name(end-2:end);
    if strcmp(postfix,'bim') || strcmp(postfix,'fix')
        name2 = [path '/' name(1:end-3) 'obj.mat'];
    else
        name2 = [path '/' name '_obj.mat'];
    end
else
    name2 = file;
end
T2 = [1  0   0;   0  1  0; 0   0   1]; %Scale
vertices= vertices*T2;
% Second, perform the initial parameterization
[sph_verts, new_name] = initParamCALD(vertices, faces, name2, confs);
%  name3 to new_name

% %Calculate the smoothed version
%    [vertices1, faces1]=smoothpatch(vertices, faces,1,5,1,1);
% 
%    % Show the mesh and smoothed mesh
%    figure, 
%     subplot(1,2,1), patch('Faces',faces,'Vertices',vertices,'FaceColor',[1 0 0],'EdgeAlpha',0);  view(3); camlight
%     subplot(1,2,2), patch('Faces',faces1,'Vertices',vertices1,'FaceColor',[0 0 1],'EdgeAlpha',0); view(3); camlight
% 

% Third, smooth the parameterization
%[vertices, faces, sph_verts, new_name] = smootheCALD(vertices, faces, sph_verts, name3, confs);

return;