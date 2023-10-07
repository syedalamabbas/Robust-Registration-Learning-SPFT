
%
%match objects together first
%

function [rvec, M] = match1_icp(fvec,atlas,max_d);

dg = [0 max_d]; meshsize = -3;
% purpose: generate surface data structure based on spharm
[X, fs] = surf_spharm(atlas,dg,meshsize);
[P, fs] = surf_spharm(fvec,dg,meshsize);

% align P to X
[P,M] = align_icp(P,X);
rvec = (M(1:3,1:3)*fvec')';
rvec(1,:) = fvec(1,:) + M(1:3,4)'*2*sqrt(pi); % Y00 = 1/(2*sqrt(pi))

rmsd(1) = SPHARM_rmsd(fvec, atlas);
rmsd(2) = SPHARM_rmsd(rvec, atlas);

if rmsd(1)< rmsd(2)
    rvec = fvec; 
end

return;