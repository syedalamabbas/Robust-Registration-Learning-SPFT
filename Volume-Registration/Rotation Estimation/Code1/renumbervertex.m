% Renumber vertices in faces such that vertex ids are contigous starting from 1

% The renumbering is such that the ordering heirarchy of vertices is preserved
% A vertex with an id lower than another vertex will have the renumbered id lower as well
% e.g. If original vertices are 30,40,60 they will be renumbered 1,2,3

% Only renumbered face connectivity is output.
% Obtain renumbered vertex crd outside the function using
% vertices_ren=vertices((orig_vertex_id),:)  where vertices=original vertex crd

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%     Copyright (C) 2010  Ravi Soni www.hermesacademy.com ravi@hermesacademy.com
%
%     Program "Voxel Surface" is a free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     Program consists of the following files
%     run					Run a demo example
%     voxel_bnd_faces       Main function
%	  voxel_vtx             
%     renumbervertex
%     plotsurf              

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [faces_ren,orig_vertex_id]=renumbervertex(faces)

% Input
% faces(nfaces,nfv)         Connectivity of faces : nfv=no. vertices of each face

% Output
% faces_ren(nfaces,nfv)     Connectivity of faces with vertices renumbered
% orig_vertex_id(nv,1)      Original vertex ids : used to obtain original vertex id from renumbered id

% Local
% nv                        No. vertices
% flag_vtx(:,1)             Flag for determining if index is a vertex id
%                           0 = index is not a vertex id
%                           1 = index is a vertex id
% ren_vertex_id(:,1)        Renumbered vertex ids : used to obtain renumbered vertex id from original id
%                           0 = index is not a vertex id
%                           ~0= renumbered id

% -------------------------------------------------------------------------

[nfaces,nfv]=size(faces);

faces_ren=[];
orig_vertex_id=[];

% Determine vertex ids in faces
for i=1:nfv
    flag_vtx(faces(:,i),1)=ones(nfaces,1);
end
flag_vtx;

orig_vertex_id=find(flag_vtx);  % non-zero entries correspond to vertex ids

nv=length(orig_vertex_id);  % no. vertices

ren_vertex_id(orig_vertex_id,1)=(1:nv);  % renumber vertices contigously starting from 1

% Determine faces with vertices renumbered
for i=1:nfv
    renvtx=ren_vertex_id(faces(:,i));
    faces_ren(:,i)=renvtx;
end
faces_ren;

% -------------------------------------------------------------------------
% End of Function