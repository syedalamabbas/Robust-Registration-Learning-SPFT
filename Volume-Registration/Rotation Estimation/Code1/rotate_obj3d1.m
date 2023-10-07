function [vertices,Tr] = rotate_obj3d1(vertices,theta_r)


T2yr = [cosd(theta_r.y)  0     sind(theta_r.y)  ;   0  1 0   ;  -sind(theta_r.y)   0   cosd(theta_r.y)  ]; %Round Y
T2zr = [cosd(theta_r.z)  -sind(theta_r.z) 0  ;sind(theta_r.z)   cosd(theta_r.z)  0 ;0  0   1  ]; %
T2xr = [1  0  0 ; 0  cosd(theta_r.x)   -sind(theta_r.x) ; 0  sind(theta_r.x)   cosd(theta_r.x) ];  % Around X
Tr=T2zr * T2yr * T2xr;    %ZYX
%  T1 = [1  0  0  ; 0 1  0  ; 0  0  1  ; -blob_center ];
%  T3 = [1  0  0  ; 0 1  0  ; 0  0  1  ; blob_center ];
% Tr=T1 * T2 * T3;
vertices=vertices * Tr;

%% (3) Display input objects (optional)
%     %Available values for Space- 'object';'param';'both'
% dispConfs.Space = 'object';
%     % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
%     %    'icosa3';'icosa4';'icosa5';'icosa6'
% dispConfs.Mesh = 'orig';  
%     % Available values for Shape- 'solid';'mesh';'both'
% dispConfs.Shade = 'both';
%     % Available values for Overlay- 'none';'adc_paramap'
% dispConfs.Overlay = 'none';
%     % Available values for Export- 'screen';'png';'both'
% dispConfs.Export = 'png';
% dispConfs.Degree = [];
% dispConfs.Template = '';
% 
% SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);
return;


