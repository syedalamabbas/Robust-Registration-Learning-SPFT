% Create surfaces from 3-D voxel data

% A voxel is considered to be a cuboid of dimensions sf(1) x sf(2) x sf(3).
% sf are the voxel scale factors for x,y,z dimensions.

% Surface of a domain corresponding to a given voxel id is comprised of the boundary faces of the voxel cuboids.
% A face is a boundary face if it belongs to only one voxel.

% Faces are triangular and faces/vertices are numbered contigously starting from 1.

% Face connectivity is such that outward normals are consistent and point outwards.
% Normal is given by a x b, where a=edge 1-2  b=edge 2-3.

% Facets are perfect with no quality issues.

% Normal consistency check and other quality checks are therefore not required.

% This function replaces MATLAB isosurface, isocap and isonormals fncs.
% MATLAB fncs do not render exact surfaces (even without smooth3).
% The surface facets have quality issues for complex shapes - bad connectivity and duplicate faces.
% Normals are not consistent and isonormals fnc runs out of memory for large data. 
% They are also cpu intensive (esp. if extra quality checks are added).

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
%     Program consists of the following files
%     run					Run a demo example
%     voxel_bnd_faces       Main function
%	  voxel_vtx             
%     renumbervertex
%     plotsurf              

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [bndfaces_ren,bndvtxcrd]=voxel_bnd_faces(voxels,sf,offset,id)

% Input
% voxels(nrow,ncol,nslice)  Ids of all voxels : nrow=ny ncol=nx nslice=nz
% sf(1,3)                   Voxel scale factors for x,y,z dimensions
%                           sf=[1 1 1] if problem does not have any dimensional scale
% offset(1,3)               Offset of crd of center of voxel(1,1,1) from global origin
%                           offset=[0 0 0] if global origin is the center of voxel(1,1,1)
% id                        Voxel id whose boundary is to be determined

% Output
% bndfaces_ren(nbndfaces,3) Connectivity of boundary faces : vertices are numbered contigously from 1 to nbndvtx
% bndvtxcrd(nbndvtx,dim)    Crd of boundary vertices : dim=3

% Local
% nrow,ncol,nslice          Range of voxel indices

% nvertices                 Total no. vertices in all voxels
% vertices(nvertices)       Vertex ids numbered 1:nvertices
% vertices_crd(nvertices,3) Crd of vertices

% bndfaces(nbndfaces,3)     Connectivity of boundary faces : vertices correspond to global ids and are not contigous 

% j,i,k                     Voxel indices 
% vtx(nv)                   Global ids of vertices of  voxel(i,j,k) : nv=8
% vtxcrd(nv,3)              Crd of vertices of  voxel(i,j,k)

% Note:
% Voxel data is ordered by row,col,slice.
% Thus indices j,i,k correspond to row,col,slice which correspond to y,x,z.
% Thus i,j,k correspond to x,y,z.

% vertices_crd only has crd of vertices for voxels that match the given id.
% If it is required to obtain all crd, move call to fnc=voxel_vtx before "if (voxels(i,j,k)==id)".

% -------------------------------------------------------------------------

[nrow,ncol,nslice]=size(voxels);

% First index=row=y   Second index=col=x
nx=ncol;
ny=nrow;
nz=nslice;

nvertices=(nx+1)*(ny+1)*(nz+1);
vertices_crd=zeros(nvertices,3);

nbndfaces=0;
bndfaces=[];

for k=1:nz
    k
    for j=1:ny
        for i=1:nx
            % Voxel(j,i,k) 
               
            % Check for boundary faces if voxel is of given id
            if (voxels(j,i,k)==id)
                % Voxel corresponds to given id 
                
                % Determine global ids and crd of voxel vertices
                [vtx,vtxcrd]=voxel_vtx(nx,ny,nz,sf,offset,i,j,k);
                
                % Store crd of vertices
                vertices_crd(vtx,:)=vtxcrd;
                
                % ---------------------------------------------------------
                % Check if x-faces are boundary faces
                if (i==1 | voxels(j,i-1,k)~=id)
                    % x-face 1 is a boundary face
                    
                    % Add tria 1-5-7 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(1) vtx(5) vtx(7)] ; 
                    bndfaces(nbndfaces,:)=face;
                    
                    % Add tria 7-3-1 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(7) vtx(3) vtx(1)];
                    bndfaces(nbndfaces,:)=face;
                end
                
                if (i==nx | voxels(j,i+1,k)~=id)
                    % x-face 2 is a boundary face
                    
                    % Add tria 2-4-8 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(2) vtx(4) vtx(8)] ; 
                    bndfaces(nbndfaces,:)=face;
                    
                    % Add tria 8-6-2 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(8) vtx(6) vtx(2)];
                    bndfaces(nbndfaces,:)=face;
                end  
                
                % ---------------------------------------------------------
                % Check if y-faces are boundary faces
                if (j==1 | voxels(j-1,i,k)~=id)
                    % y-face 1 is a boundary face
                    
                    % Add tria 1-2-6 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(1) vtx(2) vtx(6)] ; 
                    bndfaces(nbndfaces,:)=face;
                    
                    % Add tria 6-5-1 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(6) vtx(5) vtx(1)];
                    bndfaces(nbndfaces,:)=face;
                end
                
                if (j==ny | voxels(j+1,i,k)~=id)
                    % y-face 2 is a boundary face
                    
                    % Add tria 3-7-8 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(3) vtx(7) vtx(8)] ; 
                    bndfaces(nbndfaces,:)=face;
                    
                    % Add tria 8-4-3 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(8) vtx(4) vtx(3)];
                    bndfaces(nbndfaces,:)=face;
                end 
                
                % ---------------------------------------------------------
                % Check if z-faces are boundary faces
                if (k==1 | voxels(j,i,k-1)~=id)
                    % z-face 1 is a boundary face
                    
                    % Add tria 1-3-4 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(1) vtx(3) vtx(4)] ; 
                    bndfaces(nbndfaces,:)=face;
                    
                    % Add tria 4-2-1 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(4) vtx(2) vtx(1)];
                    bndfaces(nbndfaces,:)=face;
                end
                
                if (k==nz | voxels(j,i,k+1)~=id)
                    % z-face 2 is a boundary face
                    
                    % Add tria 5-6-8 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(5) vtx(6) vtx(8)] ; 
                    bndfaces(nbndfaces,:)=face;
                    
                    % Add tria 8-7-5 to boundary faces
                    nbndfaces=nbndfaces+1;
                    face=[vtx(8) vtx(7) vtx(5)];
                    bndfaces(nbndfaces,:)=face;
                end          
                
            end  % if (voxels(j,i,k)==id)
        end  % for i=1:nx
    end  % for j=1:ny
end  % for k=1:nz

clear voxels
bndfaces;
nbndfaces;
size_bndfaces=size(bndfaces);

% Renumber boundary faces such that vertex id are contigous starting from 1
[bndfaces_ren,orig_vertex_id]=renumbervertex(bndfaces);

% Obtain crd of boundary vertices
vertices=1:nvertices;
bndvtx=vertices(orig_vertex_id);
bndvtxcrd=vertices_crd(bndvtx,:);
size_bndvtxcrd=size(bndvtxcrd);
                
% -------------------------------------------------------------------------
% End of Function