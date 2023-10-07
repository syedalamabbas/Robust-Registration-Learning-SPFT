% ============================================
%
% Goal: decompistion P to X 
%


function [R,qT] = decomp_eigen(P,Xc)


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
   delta = [A(2,3) ; A(3,1) ; A(1,2)];
   Q     = [trace(sig)  delta' ; delta  sig+sig'-trace(sig)*eye(3)];
   [V,D] = eig(Q);
   [val,ind] = max(diag(D));
   qR    = V(:,ind); % optimal rotation
   q0    = qR(1);
   q1    = qR(2);
   q2    = qR(3);
   q3    = qR(4);
   R     = [q0^2+q1^2-q2^2-q3^2  2*(q1*q2 - q0*q3)  2*(q1*q3 + q0*q2) ; ...
	    2*(q1*q2 + q0*q3)  q0^2+q2^2-q1^2-q3^2  2*(q2*q3 - q0*q1) ; ...
	    2*(q1*q3 - q0*q2) 2*(q2*q3 + q0*q1)  q0^2+q3^2-q1^2-q2^2];
   qT    = muXc' - R*muP'; % optimal translation
   % fprintf( 'R=' );
%     R
%     [thetaX, thetaY, thetaZ] = factor_rot_xyz(R');
%       disp(sprintf('Factor rotation [theta]xyz: %0.2f %0.2f %0.2f',round(thetaX/pi*180),round(thetaY/pi*180),round(thetaZ/pi*180)+0));
