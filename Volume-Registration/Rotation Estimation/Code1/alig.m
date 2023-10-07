

function [rmsd]=alig(filename, confs)

% global fact;
% global sgm;
% gran = confs.GroupAlpha;  % gran = # of alphas to process together
% res = [confs.BaseRes confs.HierarchyStep confs.HierarchyDepth confs.Top_K confs.GammaRes]; 
% base res R + step of hierarchy Hs + depth of hierarchy Hd + top N + 3rd Angle res (gammares)


% Load a template object
load(confs.Template);
tfvec = fvec;

load(filename);

max_d=confs.MaxSPHARMDegree;
% rmsd_org = SPHARM_rmsd(fvec, tfvec);
% 
% if rmsd_org > 0

    % create samples in rotation space
%     [alpha, beta, gamma] = utl_eas(res); % euler_angle hierarchy
% % factorial(170) = Inf
%        for i=0:170 
%             fact(i+1) = factorial(i);
%        end
%        utl_sgm(15);
    % initial alignment using ICP
%     [fvec_icp, Robj_icp] = match1_icp(fvec, tfvec, max_d);
%     [fvec_icp, Aprm] = match_param_hie(fvec_icp,tfvec,alpha,beta,gamma,gran,res,max_d);
%     [fvec_icp, Robj_icp,rmsd] = match_cps(fvec_icp, tfvec, max_d);

    % assume that fvec and atlas are roughly aligned in the object space
    [rmsd] = match_cps(fvec, tfvec, max_d);

    %rmsd_icp = SPHARM_rmsd(fvec_icp, tfvec);
    %rmsd_cps = SPHARM_rmsd(fvec_cps, tfvec);

%     disp(sprintf('RMSD: org %0.3f, icp %0.3f, cps %0.3f',rmsd_org,rmsd,rmsd_cps))

%     if rmsd_org > min(rmsd_icp,rmsd_cps)    
%         if rmsd_icp < rmsd_cps
% %             fvec = fvec_icp;
% %             Robj = Robj_icp;
%               rmsd=rmsd_icp;
%             
%         else
% %             fvec = fvec_cps;
% %             Robj = Robj_cps;
%               rmsd=rmsd_cps;
%         end
%     end
%     
    % assume that fvec and atlas are roughly aligned in both object and
    % parameter spaces
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
%         rmsd = SPHARM_rmsd(fvec, tfvec);
% 
%        % disp(sprintf('Base Res %d, Step %d, Depth %d, Top %d, Gamma %d: RMSD %0.3f => %0.3f => %0.3f',res,rmsd));
% 
%         if sum(abs(Aprm(:)))==0
%             break;
%         end
%     end
% else 
%     rmsd=0;
%   %  rmsd_cps=0;
%end


return;
