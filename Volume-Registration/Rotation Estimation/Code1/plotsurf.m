% Plot faceted surface

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

function [fighndl]=plotsurf(faces,vertex,facecolor,edgecolor,name,tit,fighndl)

% Input
% faces(nf,3)           Face connectivity of isosurface
% vertex(nv,3)          Vertex coordinates of isosurface
% facecolor(3)          RGB value for faces
% edgecolor(3)          RGB value for edges
% name                  Name of image
% fighndl               Figure handle of existing figure

% Output
% fighndl               Figure handle 

% Local

% -------------------------------------------------------------------------

value=exist('fighndl','var');

if (value==0 | (isempty(fighndl) | fighndl>=0))
    % Plot surfaces if fighndl is undefined or is non-negative
    if (~isempty(faces))
        % Display Isosurface
        if (value==0 | isempty(fighndl))
            fighndl=figure('Name',name);  % new figure 
            title(tit);
        else
            figure(fighndl);  % existing figure
        end

        hiso = patch('Faces',faces,'Vertices',vertex,... 
                     'FaceColor',facecolor,'EdgeColor',edgecolor);            

        % lightangle(305,30); % from MATLAB example
        lightangle(305,15);  % gives better results
        set(gcf,'Renderer','zbuffer'); 
        lighting phong
        view(45,30) 
        axis tight 
        daspect([1,1,1])
        grid;
    else
        % No isosurface : do not plot
        fighndl=[];
    end
end

% -------------------------------------------------------------------------
% End of Function