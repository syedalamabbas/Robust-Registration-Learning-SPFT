function [theta] = factor_rot_xyz(R)
R=real(R);
thetaY = asin(R(1,3));
if (thetaY < pi/2)
    if (thetaY > -pi/2)
        thetaX = atan2(-R(2,3),R(3,3));  % Rotation around z
        thetaZ = atan2(-R(1,2),R(1,1));
    else
        disp('WARNING (Factor Rotation): thetaY = -pi/2,  set thetaZ = 0');
        thetaX = -atan2(R(2,1),R(2,2));
        thetaZ = 0;
    end
else
    disp('WARNING (Factor Rotation): thetaY = pi/2, set thetaZ = 0');
    thetaX = atan2(R(2,1),R(2,2));
    thetaZ = 0;
end

    theta.x=abs(round(thetaX/pi*180));
    theta.y=abs(round(thetaY/pi*180));
    theta.z=abs(round(thetaZ/pi*180));
return;