
function [rmsd]=eigenvector(filename, confs)

 scale = 0;

% Load a template object
load(confs.Template);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale); % fix the length of coefficent template
tfvec = fvec;

load(filename);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);% fix the length of coefficent test object

if ~exist('faces', 'var') || ~exist('vertices', 'var') || ~exist('sph_verts', 'var') || ~exist('fvec', 'var')
    disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
    return;
end

[rmsd(1,:)] = match(fvec, tfvec, max_d);  %icp
 [rmsd(1,:)] = match_Eigenvector(fvec,tfvec,max_d);
  %[rmsd] = SPHARM_rmsd(fvec, tfvec,max_d);

return;
