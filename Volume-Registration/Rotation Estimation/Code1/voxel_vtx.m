% Determine global ids of vertices of a given voxel and their crd

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

function [vtx,vtxcrd]=voxel_vtx(nx,ny,nz,sf,origin,i,j,k)

% Input
% nx,ny,nz          Range of voxel indices along x,y,z
% sf(1,3)           Voxel scale factors for x,y,z dimensions
% origin(1,3)       Origin : offset of crd of voxel(1,1,1) from global origin
% i,j,k             Voxel indices along x,y,z

% Output
% vtx(nv)           Global ids of voxel vertices : nv=8
% vtxcrd(nv,3)      Crd of vertices

% Local

% Note:
% Indices i,j,k correspond to x,y,z.
% However voxel data is ordered by row,col,slice.
% Thus i corresponds to col, j to row, z to slice.
% This should be accounted for in call to function.

% -------------------------------------------------------------------------

% Global ids of voxel vertices
vtx(1)=(nx+1)*(ny+1)*(k-1) + (nx+1)*(j-1) + i;
vtx(2)=vtx(1)+1;
vtx(3)=(nx+1)*(ny+1)*(k-1) + (nx+1)*(j) + i;
vtx(4)=vtx(3)+1;
vtx(5)=(nx+1)*(ny+1)*(k) + (nx+1)*(j-1) + i;
vtx(6)=vtx(5)+1;
vtx(7)=(nx+1)*(ny+1)*(k) + (nx+1)*(j) + i;
vtx(8)=vtx(7)+1;

% Crd of center of voxel
indx=[i-1 j-1 k-1];
rc=origin+sf.*indx+sf/2;

% Crd of voxel vertices
knt=0;
for kk=1:2
    for jj=1:2
        for ii=1:2
            knt=knt+1;
            delta=[(-1)^ii*sf(1) (-1)^jj*sf(2) (-1)^kk*sf(3)]/2;
            vtxcrd(knt,:)=rc+delta;
        end
    end
end
vtxcrd;

% -------------------------------------------------------------------------
% End of Function