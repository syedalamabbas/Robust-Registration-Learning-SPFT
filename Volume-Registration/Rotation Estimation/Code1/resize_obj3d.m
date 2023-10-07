
function [vertices] = resize_obj3d(vertices,scale)


TScale = [scale  0     0  ;   0  scale 0   ;  0  0   scale  ]; %Scale
vertices=vertices * TScale;

return;


