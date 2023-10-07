function [thetaX1, thetaY1, thetaZ1] = factor_rot_xyz2(R)

if ((R(1,3) ~= 1)&& (R(1,3) ~= -1))
 
thetaY1 = -asin(R(1,3));
thetaY2 = pi - thetaY1;

% if (thetaY < pi/2)
%     if (thetaY > -pi/2)
  thetaX1 = atan2((R(2,3)/cos(thetaY1)),(R(3,3)/cos(thetaY1)));  % Rotation around z
  thetaX2 = atan2((R(2,3)/cos(thetaY2)),(R(3,3)/cos(thetaY2)));
  thetaZ1 = atan2((R(1,2)/cos(thetaY1)),(R(1,1)/cos(thetaY1)));
  thetaZ2 = atan2((R(1,2)/cos(thetaY2)),(R(1,1)/cos(thetaY2)));
  
else
     thetaZ1 = 0;
     
      if (R(1,3) == -1)
         thetaY1=pi/2;
         thetaX1 = thetaZ1+atan2(R(2,1),R(3,1));
     else
       thetaY1=-pi/2;
       thetaX1 = -thetaZ1+ atan2(-R(2,1),-R(3,1));
     end

end
disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',(thetaX1/pi*180),(thetaY1/pi*180),(thetaZ1/pi*180)));
disp(sprintf('Factor rotation  [theta]xyz: %0.2f %0.2f %0.2f',(thetaX2/pi*180),(thetaY2/pi*180),(thetaZ2/pi*180)));

return;

