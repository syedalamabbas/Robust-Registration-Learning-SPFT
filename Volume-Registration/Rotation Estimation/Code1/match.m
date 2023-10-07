
%
%match objects together first
%

function [rmsd,theta] = match(fvec,atlas,max_d,confs)

dg = [0 max_d]; 
meshsize = -3;
% purpose: generate surface coefficents data structure based on spharm
 [X, fs] = surf_spharm(atlas,dg,meshsize);
 [P, fs] = surf_spharm(fvec,dg,meshsize);
%P = real(fvec);
%X = real(atlas);

[P,M] = match_eigen(P,X);
%[P,M] = match_eigen((fvec),(atlas));
rvec = (M(1:3,1:3)*fvec')';
%rvec(1,:) = fvec(1,:) + M(1:3,4)'*2*sqrt(pi); % Y00 = 1/(2*sqrt(pi))
 R=M(1:3,1:3);
 [thetaX, thetaY, thetaZ] = factor_rot_xyz(R');
 %disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',abs(round(thetaX/pi*180)),abs(round(thetaY/pi*180)),abs(round(thetaZ/pi*180))));
theta.x=abs(round(thetaX/pi*180));
theta.y=abs(round(thetaY/pi*180));
theta.z=abs(round(thetaZ/pi*180));
% TL=M(1:3,4);
rmsd1 = SPHARM_rmsd(fvec, atlas);
rmsd2 = SPHARM_rmsd(rvec, atlas);
rmsd=rmsd2;
% if( rmsd1<rmsd2)
%     rmsd= rmsd1; 
% end

return;