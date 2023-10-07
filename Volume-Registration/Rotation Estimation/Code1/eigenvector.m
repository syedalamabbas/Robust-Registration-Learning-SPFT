
function [rmsd,theta, R]=eigenvector(filename, confs)

 scale = 0;

% Load a template object
load(confs.Template);

%TL=length(vertices(:,1))
 [fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale); % fix the length of coefficent template

tfvec = fvec;
%tfvec= real(fvec(2:4,:));  



load(filename);
%FL=length(vertices(:,1))
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);% fix the length of coefficent test object
% ffvec=real(fvec(2:4,:));
% max_d=1;

if ~exist('faces', 'var') || ~exist('vertices', 'var') || ~exist('sph_verts', 'var') || ~exist('fvec', 'var')
    disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
    return;
end

%[rmsd(1,:), theta] = match(fvec, tfvec, max_d,confs);  %icp
 [rmsd(1,:),theta, R] = match_Eigenvector(fvec,tfvec,max_d,confs);
  %[rvec, M] = match_cpsO(fvec,tfvec,max_d);
%  tvertives=vertives; 
%  Rvertices=(real(R)*vertices')';
%  [Rvertices, max_d] = fixed_fvec(Rvertices,confs.MaxSPHARMDegree,scale); % fix the length of coefficent template
%  [tvertives, max_d] = fixed_fvec(tvertives,confs.MaxSPHARMDegree,scale); % fix the length of coefficent template
%  
%  [rmsdO(1,:)] = SPHARM_rmsd(Rvertices, tvertives);
%[rmsdO(1,:)]=0;
return;
