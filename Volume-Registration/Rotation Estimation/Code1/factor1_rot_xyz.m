% INVRPY - inverse of Roll Pitch Yaw transform
%
% Usage:  [rpy1, rpy2] = invrpy(RPY)
% 
% Argument:  RPY - 4x4 Homogeneous transformation matrix or 3x3 rotation matrix
% Returns:  rpy1 = [phi1, theta1, psi1] - the 1st solution and
%           rpy2 = [phi2, theta2, psi2] - the 2nd solution
%
%  rotx(phi1)*roty(theta1)*rotz(psi1) = RPY
%


function [rpy1, rpy2] = factor1_rot_xyz(RPY)

    phi1 = atan2(RPY(2,1), RPY(1,1));
    phi2 = phi1 + pi;
    
    theta1 = atan2(-RPY(3,1), cos(phi1)*RPY(1,1) + sin(phi1)*RPY(2,1));
    theta2 = atan2(-RPY(3,1), cos(phi2)*RPY(1,1) + sin(phi2)*RPY(2,1));
    
    psi1 = atan2(sin(phi1)*RPY(1,3) - cos(phi1)*RPY(2,3), ...
                 -sin(phi1)*RPY(1,2) + cos(phi1)*RPY(2,2));
    psi2 = atan2(sin(phi2)*RPY(1,3) - cos(phi2)*RPY(2,3), ...
                 -sin(phi2)*RPY(1,2) + cos(phi2)*RPY(2,2));
    
    rpy1 = [phi1/pi*180, theta1/pi*180, psi1/pi*180];
    rpy2 = [phi2/pi*180, theta2/pi*180, psi2/pi*180];