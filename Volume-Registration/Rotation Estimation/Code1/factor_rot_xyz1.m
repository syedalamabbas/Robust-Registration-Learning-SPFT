function [thetaX, thetaY, thetaZ] = factor_rot_xyz1(R)

thetaY = asin(R(1,3));
if (thetaY < pi/2)
    if (thetaY > -pi/2)
        thetaX = -atan2(-R(2,3),R(3,3));  % Rotation around z
        thetaZ = -atan2(-R(1,2),R(1,1));
       else
        disp('WARNING (Factor Rotation): thetaY = -pi/2,  set thetaZ = 0');
        thetaX = atan2(R(2,1),R(2,2));
        thetaZ = 0;
    end
else
    disp('WARNING (Factor Rotation): thetaY = pi/2, set thetaZ = 0');
    thetaX = -atan2(R(2,1),R(2,2));
    thetaZ = 0;
end
 disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',(thetaX/pi*180),(thetaY/pi*180),(thetaZ/pi*180)));
% if (thetaY < 0)
%     thetaY=-thetaY;
%     thetaX=-(thetaX+pi);
%     thetaZ= -(pi+thetaZ);
% end
    

return;

