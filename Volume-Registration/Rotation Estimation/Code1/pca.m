

function [rmsd1,theta, R]=pca(filename, confs)

load(confs.Template);
%[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,0);
tvertices = vertices; %tsph_verts=sph_verts; tfvec = fvec;
tfvec = fvec;

load(filename);
% [fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,0);
[path,name,ext] = fileparts(filename);
% 
% if ~exist('faces', 'var') | ~exist('vertices', 'var') | ~exist('sph_verts', 'var') | ~exist('fvec', 'var')
%     disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
%     return;
% end

% center P and X first
[P,cmP] = center(tvertices);
[X,cmX] = center(vertices);

Xc = X;

muP   = mean(P);
muXc  = mean(Xc);   
N     = length(P);
sig   = zeros(3);
for k = 1 : N
  sig = sig + (P(k,:)' * Xc(k,:));
end
sig   = sig/N;
sig   = sig - (muP' * muXc); % cross-covariance matrix
A     = sig - sig';
% delta = [A(2,3) ; A(3,1) ; A(1,2)];
% Q     = [trace(sig)  delta' ; delta  sig+sig'-trace(sig)*eye(3)];
% [V,D] = eig(Q);
% [val,ind] = max(diag(D));
% R    = V(:,ind); % optimal rotation

[U,S,V] = svd(sig);
R=V;
if (det(R)<0)
    disp('WARNING (Parameter Space): rotoinversion!! Change back to pure rotation');
    R(2,:) = R(2,:)*(-1);
end
rvec = tfvec*R; 
%rvertices = tvertices*R'; 

rmsd1(1,:) = SPHARM_rmsd(rvec, fvec);



[theta] = factor_rot_xyz(R);


return;
