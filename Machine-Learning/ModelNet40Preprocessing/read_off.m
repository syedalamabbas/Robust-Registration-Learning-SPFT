function [vertex,face] = read_off(filename)

% read_off - read data from OFF file.
%
% [vertex,face] = read_off(filename);
%
% 'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
% 'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
% Copyright (c) 2003 Gabriel Peyr�

fid = fopen(filename,'r');
if( fid==-1 ) 
    error('Can''t open the file.'); 
    return;
end

str = fgets(fid); % -1 if eof
if ~strcmp(str(1:3), 'OFF') 
    error('The file is not a valid OFF one.'); 
end

% Added support for same line of nvertices and nfaces, Syed Alam Abbas
% 9/17/2022
if(length(str) ~= 4)
    str = str(4:length(str));
else
    str = fgets(fid);
end

[a,str] = strtok(str); 
nvert = str2double(a);
[a,~] = strtok(str); 
nface = str2double(a);


[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert 
    warning('Problem in reading vertices.');
end

A = reshape(A, 3, cnt/3);
vertex = A;

% read Face 1 1088 480 1022
[A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
if cnt~=4*nface 
    warning('Problem in reading faces.');
end

A = reshape(A, 4, cnt/4);
face = A(2:4,:)+1;

fclose(fid);