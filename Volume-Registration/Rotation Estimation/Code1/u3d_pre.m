function [fvc3,fvc4] = u3d_pre
% Syntax: 
%     fvc = u3d_pre
%
% Description:
%    This function generates the input for the MESH_TO_LATEX function from
%    Alexandre Gramfort from your surface-graphs. The surface graphs 3d-model can be
%    displayed in u3d - format in pdf or if you don't delete the u3d-files
%    which are generated by MESH_TO_LATEX you can use deepview from
%    righthemisphere(http://www.righthemisphere.com/products/client-products/deep-view) 
%    to embed the modell in Microsoft-Office products (Word, Excel, PowerPoint). 
%
%    The output is a structure generated by the standard MATLAB function
%    surf2patch which can used in MESH_TO_LATEX: 
%       fvc.vertices         -->    points
%       fvc.faces            -->    faces
%       fvc.facevertexcdata  -->    face_vertex_data
%
% Example:
%       [X,Y,Z] = peaks(30);
%       surf(X,Y,Z); 
%       fvc = u3d_pre; 
%       mesh_to_latex('tex/mesh',fvc.vertices,uint32(fvc.faces),fvc.facevertexcdata);
%
% Bugs and suggestions:
%    Please send to Sven Koerner: koerner(underline)sven(add)gmx.de
%
% See also:
%    MESH_TO_LATEX on the File Exchange
%    http://www.righthemisphere.com/products/client-products/deep-view
% 
%
% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.
%
% Programmed by Sven Koerner: koerner(underline)sven(add)gmx.de
% Date: 2010/04/14



% search for surface-graph 
h = findobj('type','surface');

if ~isempty(h)    
    % get defined data-points
    X = get(h(end,1), 'XData');
    Y = get(h(end,1), 'YData');
    Z = get(h(end,1), 'ZData');
    
    % scaled color to unscaled r
    cdata = get(h(end,1), 'CData');
    siz = size(cdata);
    cmap = colormap;
    nColors = size(cmap,1);
    cax = caxis;
    idx = ceil( (double(cdata) - cax(1)) / (cax(2)-cax(1)) * nColors);
    idx(idx<1) = 1;
    idx(idx>nColors) = nColors;
    %handle nans in idx
    nanmask = isnan(idx);
    idx(nanmask)=1; %temporarily replace w/ a valid colormap index
    realcolor = zeros(siz);
    for i = 1:3,
        c = cmap(idx,i);
        c = reshape(c,siz);
        realcolor(:,:,i) = c;
    end
    
    fvc3 = surf2patch(X,Y,Z,realcolor, 'triangles' );
    fvc4 = surf2patch(X,Y,Z,realcolor);
else
    fvc3.vertices            = [];
    fvc3.faces               = [];
    fvc3.facevertexcdata     = [];
    disp('No surface graph found.');
     fvc4.vertices            = [];
    fvc4.faces               = [];
    fvc4.facevertexcdata     = [];
    disp('No surface graph found.');

end;
