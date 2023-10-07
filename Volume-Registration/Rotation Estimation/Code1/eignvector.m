
function [fvec,confs]=eignvector(filename, confs)

scale = 0;

% Load a template object
load(confs.Template);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);
%tvertices = vertices; tsph_verts=sph_verts; tfvec = fvec;
tfvec = fvec;

load(filename);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);

if ~exist('faces', 'var') || ~exist('vertices', 'var') || ~exist('sph_verts', 'var') || ~exist('fvec', 'var')
    disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
    return;
end

rmsd_org = SPHARM_rmsd(fvec, tfvec);

if rmsd_org > 0

    [fvec_icp, Robj_icp] = match_icp(fvec, tfvec, max_d);

     [theta1X, theta1Y, theta1Z] = factor_rot_xyz(Robj_icp');
      disp(sprintf('Factor rotation 1 for %s  [theta]xyz: %0.2f %0.2f %0.2f',filename,round(theta1X/pi*180),round(theta1Y/pi*180),round(theta1Z/pi*180)+0));
end


return;
