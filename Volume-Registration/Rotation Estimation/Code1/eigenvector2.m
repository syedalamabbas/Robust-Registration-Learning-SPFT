

%function [vertices, sph_verts, faces, fvec, new_name]=align_CPS_SHREC2(filename, confs)
function [rmsd1,theta, R]=eigenvector2(filename, confs)
% global fact;
% global sgm;
% 
% gran = confs.GroupAlpha;  % gran = # of alphas to process together
% res = [confs.BaseRes confs.HierarchyStep confs.HierarchyDepth confs.Top_K confs.GammaRes]; 
% % base res R + step of hierarchy Hs + depth of hierarchy Hd + top N + 3rd Angle res (gammares)
% 
% switch upper(deblank(char(confs.NormalizeSize)))
%     case 'YES'
%         scale = 1;
%     case 'NO'
        scale = 0;
% end
% 
% % factorial(170) = Inf
% for i=0:170 
%     fact(i+1) = factorial(i);
% end
% 
% utl_sgm(15);

% Load a template object
load(confs.Template);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);
%tvertices = vertices; tsph_verts=sph_verts; tfvec = fvec;
tfvec = fvec;

load(filename);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);
[path,name,ext] = fileparts(filename);

if ~exist('faces', 'var') | ~exist('vertices', 'var') | ~exist('sph_verts', 'var') | ~exist('fvec', 'var')
    disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
    return;
end

rmsd_org = SPHARM_rmsd(fvec, tfvec);
rmsd1(1,:)=rmsd_org;
% if rmsd_org > 0
% 
%     % create samples in rotation space
%     [alpha, beta, gamma] = utl_eas(res); % euler_angle hierarchy
% 
%     % initial alignment using ICP
%     [fvec_icp, Robj_icp] = match_icp(fvec, tfvec, max_d);
%     [fvec_icp, Aprm] = match_param_hie(fvec_icp,tfvec,alpha,beta,gamma,gran,res,max_d);
%     [fvec_icp, Robj_icp] = match_cps(fvec_icp, tfvec, max_d);
%     
%    % [fvec_icp, Robj_icp] = match_cpsT(fvec_icp, tfvec, max_d);
%     % assume that fvec and atlas are roughly aligned in the object space
    [fvec_cps, Robj_cps] = match_cps(fvec, tfvec, max_d);
  RR=Robj_cps;
%     rmsd_icp = SPHARM_rmsd(fvec_icp, tfvec);
    rmsd_cps = SPHARM_rmsd(fvec_cps, tfvec);

 R=Robj_cps(1:3,1:3);
 [theta] = factor_rot_xyz(R');


%     disp(sprintf('RMSD: org %0.3f, icp %0.3f, cps %0.3f',rmsd_org,rmsd_icp,rmsd_cps))
% 
%     if rmsd_org > min(rmsd_icp,rmsd_cps)    
%         if rmsd_icp < rmsd_cps
%             fvec = fvec_icp;
%             Robj = Robj_icp;
%         else
%             fvec = fvec_cps;
%             Robj = Robj_cps;
%         end
%     end
%     
%     % assume that fvec and atlas are roughly aligned in both object and
%     % parameter spaces
%     for k = 1:5
%         rmsd(1) = SPHARM_rmsd(fvec, tfvec);
% 
%         % use cps to align objects together first
%         [fvec_a, Robj] = match_cps(fvec, tfvec, max_d);
%         rmsd(2) = SPHARM_rmsd(fvec_a, tfvec);
% 
%         if rmsd(2)<rmsd(1)
%             fvec = fvec_a; % update
%         end    
% 
%         % assume aligned in object space, register the parameterization
%         [fvec, Aprm] = match_param_hie(fvec,tfvec,alpha,beta,gamma,gran,res,max_d);
%         rmsd(3) = SPHARM_rmsd(fvec, tfvec);
% 
%         disp(sprintf('Base Res %d, Step %d, Depth %d, Top %d, Gamma %d: RMSD %0.3f => %0.3f => %0.3f',res,rmsd));
% 
%         if sum(abs(Aprm(:)))==0
%             break;
%         end
%     end
% else
%     disp(sprintf('RMSD = %d: Individual (%s) is the same as template (%s)',rmsd_org,filename,confs.Template));
% end
% 
% new_name = sprintf('%s/%sSHREC_reg.mat',confs.OutDirectory, name);
% % if exist(new_name,'file')
% %     prompt = {'Enter new filename:'};
% %     dlg_title = 'New File Name';
% %     num_lines = 1;
% %     def = {new_name};
% %     answer = inputdlg(prompt,dlg_title,num_lines,def);    
% %     new_name = answer{1};
% % end
% save(new_name, 'fvec');

clear('fact','sgm');

return;
