


function [rmsd]= match_c(fvec,atlas,max_d)

% dg = [0 max_d]; meshsize = -3;
%generate surface data structure based on spharm
% [X, fs] = surf_spharm(atlas,dg,meshsize);
% [P, fs] = surf_spharm(fvec,dg,meshsize);

% align P to X
%[P,M] = align_cps(P,X);


[P,M] = match_eign((fvec),(atlas));
 R=M(1:3,1:3);
 [thetaX, thetaY, thetaZ] = factor_rot_xyz(R');
 if (abs(thetaX) > 1.57)
    disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',abs(round(thetaX/pi*180))-90,abs(round(thetaY/pi*180)),abs(round(thetaZ/pi*180))-90));
     disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',abs(round(thetaX/pi*180)),abs(round(thetaY/pi*180)),abs(round(thetaZ/pi*180))));
 else
 disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',abs(round(thetaX/pi*180)),abs(round(thetaY/pi*180)),abs(round(thetaZ/pi*180))));
 end
%  [euler1, euler2] = factor_rot_zyz(R');
%  euler1
%  euler2
% % disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',round(thetaX/pi*180),round(thetaY/pi*180),round(thetaZ/pi*180)));
%  [euler1, euler2] = factor1_rot_xyz(R');
%  euler1
%  euler2
rvec = (M(1:3,1:3)*fvec')';
%rvec(1,:) = fvec(1,:) + M(1:3,4)'*2*sqrt(pi); % Y00 = 1/(2*sqrt(pi))

rmsd1 = SPHARM_rmsd(fvec, atlas);
rmsd = SPHARM_rmsd(rvec, atlas);

% if rmsd1<rmsd2
%     %rvec = fvec;
%     rmsd=rmsd1;
% end

return;